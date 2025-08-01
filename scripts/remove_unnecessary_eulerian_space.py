import pyvista 
import meshio
import os, sys, glob
import subprocess
import math 
from natsort import natsorted
import multiprocessing 
import pdb
import numpy as np 
import re


def write_pvd(basename, dt, nsteps, extension, nprocs_sim=1):

    prefix = '''<?xml version="1.0"?>
    <VTKFile type="Collection" version="0.1"
             byte_order="LittleEndian"
             compressor="vtkZLibDataCompressor">
      <Collection>
    '''

    suffix = '''  </Collection>
    </VTKFile>
    '''

    initialized = False

    for n in range(nsteps):
        for proc in range(nprocs_sim):

            if not initialized:

                filename_out = basename + '.pvd'
                print("filename_out = ", filename_out)

                f_write = open(filename_out, 'w')
                f_write.write(prefix)
                initialized = True

            tmp_str = '    <DataSet timestep="'
            tmp_str += '{:.14f}'.format(dt * n)
            tmp_str += '" group="" part="'
            tmp_str += str(proc) + '"'

            tmp_str += ' file="'
            if nprocs_sim > 1:
                tmp_str += basename + str(n).zfill(4) + '/' # sorted into directories 
            tmp_str += basename + str(n).zfill(4) + '.' 
            if nprocs_sim > 1:
                tmp_str += str(proc) + '.'        
            tmp_str += extension
            tmp_str += '"/>\n'

            f_write.write(tmp_str)

    f_write.write(suffix)
    f_write.close()


def read_distributed_vtr(dir_name):

    files = natsorted(glob.glob(dir_name + "/*.vtr"))
    blocks = pyvista.MultiBlock([pyvista.RectilinearGrid(f) for f in files])
    return blocks.combine(merge_points=True, tolerance=0.0001)


def sort_points_and_point_data(points, point_data_U=None, point_data_P=None):
    '''
    Sorts points by x,y,z
    Applies permutation to point and point_data
    '''

    indices = np.lexsort((points[:,0], points[:,1], points[:,2])) 

    points_sorted = points[indices]

    if point_data_U is not None:
        point_data_sorted_U = point_data_U[indices]
    else:
        point_data_sorted_U = None

    if point_data_P is not None:
        point_data_sorted_P = point_data_P[indices]
    else:
        point_data_sorted_P = None

    return points_sorted, point_data_sorted_U, point_data_sorted_P


def generate_cells(NX,NY,NZ):
    '''
    generates cells in meshio format 
    '''

    get_1d_idx = lambda i,j,k : i + j*NX + k*NX*NY

    cells = []

    # every point gets a cell 
    # no cells on the top row of points 
    for i in range(NX-1):
        for j in range(NY-1):
            for k in range(NZ-1):

                # vtk order points the normal of the "bottom" face up 
                # and the top face also up 
                cell_tmp = []
                cell_tmp.append(get_1d_idx(i  , j  , k  ))
                cell_tmp.append(get_1d_idx(i+1, j  , k  ))
                cell_tmp.append(get_1d_idx(i+1, j+1, k  ))
                cell_tmp.append(get_1d_idx(i  , j+1, k  ))
                cell_tmp.append(get_1d_idx(i  , j  , k+1))
                cell_tmp.append(get_1d_idx(i+1, j  , k+1))
                cell_tmp.append(get_1d_idx(i+1, j+1, k+1))
                cell_tmp.append(get_1d_idx(i  , j+1, k+1))

                cells.append(cell_tmp)

    return cells 



def convert_mesh_to_center_points(mesh, NX, NY, NZ):
    ''' 
    Takes cell data mesh and converts cell_array 'U' to point_array at center 
    '''

    mesh_cell_centers = mesh.cell_centers()
    points_sorted, point_data_sorted_U, point_data_sorted_P = sort_points_and_point_data(mesh_cell_centers.points, mesh_cell_centers.point_data['U'], mesh_cell_centers.point_data['P'])
    cells_sorted = generate_cells(NX,NY,NZ)

    mesh_points_meshio_format = meshio.Mesh(points_sorted, 
                                        cells={"hexahedron": cells_sorted}, 
                                        point_data={'U': point_data_sorted_U, 'P': point_data_sorted_P})

    return pyvista.wrap(mesh_points_meshio_format)



def remove_eulerian_space(basename, 
                          boundary_mesh, 
                          nsteps, 
                          label='_restricted_cells', 
                          extension='vtu', 
                          convert_to_point_data=True, 
                          NX=None, 
                          NY=None, 
                          NZ=None, 
                          proc_num=0, 
                          nprocs=1):

    for i in range(nsteps):
        if (i % nprocs) == proc_num:

            dir_name = basename + str(i).zfill(4)

            fname_out = basename + label + str(i).zfill(4) + '.' + extension

            # read distributed vtr 
            mesh = read_distributed_vtr(dir_name)

            # if there is cell data convert it 
            if convert_to_point_data:

                if ('U' in mesh.cell_data) and ('P' in mesh.cell_data):

                    if (NX is None) or (NY is None) or (NZ is None):
                        raise ValueError("Must provide values for NX,NY,NZ when converting cell to point data")

                    mesh = convert_mesh_to_center_points(mesh, NX, NY, NZ)

                # check that there is point data if there was no cell data to convert 
                elif not (('U' in mesh.point_data) and ('P' in mesh.point_data)):
                    raise ValueError('Could not find U and P in cell or point data, point data requested')


            selected = mesh.select_enclosed_points(boundary_mesh, tolerance=1.0e-10, inside_out=False, check_surface=True)

            # remove the exterior 
            # all_scalars=True keeps cells that intersect boundary 
            # all_scalars=False keeps cells that have at least one interior point 
            mesh_inside = selected.threshold(0.5, scalars="SelectedPoints", all_scalars=False) 
            #mesh_inside = selected.threshold(0.5, scalars="SelectedPoints", all_scalars=True) 

            mesh_inside.save(fname_out)


def remove_eulerian_space_single_frame(basename, 
                          boundary_mesh, 
                          frame_number, 
                          label='_restricted_cells', 
                          extension='vtu', 
                          point_data=False, 
                          NX=None, 
                          NY=None, 
                          NZ=None):

    i = frame_number

    dir_name = basename + str(i).zfill(4)

    fname_out = basename + label + str(i).zfill(4) + '.' + extension

    # read distributed vtr 
    mesh = read_distributed_vtr(dir_name)

    fname_out_cells = basename + "no_restrict_cells" + str(i).zfill(4) + '.' + extension
    mesh.save(fname_out_cells)

    if point_data:

        if (NX is None) or (NY is None) or (NZ is None):
            raise ValueError("Must provide values for NX,NY,NZ when calling points")

        mesh_point_data = convert_mesh_to_center_points(mesh, NX, NY, NZ)
        selected = mesh_point_data.select_enclosed_points(boundary_mesh, tolerance=1.0e-10, inside_out=False, check_surface=True)

    else:
        selected = mesh.select_enclosed_points(boundary_mesh, tolerance=1.0e-10, inside_out=False, check_surface=True)

    # remove the exterior 
    # all_scalars=True keeps cells that intersect boundary 
    # all_scalars=False keeps cells that have at least one interior point 
    mesh_inside = selected.threshold(0.5, scalars="SelectedPoints", all_scalars=False) 
    #mesh_inside = selected.threshold(0.5, scalars="SelectedPoints", all_scalars=True) 

    mesh_inside.save(fname_out)


def find_box_size():

    box_size_found = False 
    NX = None 
    NY = None 
    NZ = None

    log_file = None
    if os.path.isfile('IB3d.log'):
        print("log file found, trying to open")
        log_file = open('IB3d.log', 'r')
    elif os.path.isfile('../IB3d.log'):
        log_file = open('../IB3d.log', 'r')

    if not log_file:
        print("Log file not found")
        return box_size_found, NX, NY, NZ

    log_file_string = log_file.read()

    match_x = re.search(r"(NX\s+= )(\d+)", log_file_string)
    match_y = re.search(r"(NY\s+= )(\d+)", log_file_string)
    match_z = re.search(r"(NZ\s+= )(\d+)", log_file_string)

    if match_x and match_y and match_z:

        # print("match_x = ", match_x, match_x.group(), match_x.group(2))
        # print("match_y = ", match_y, match_y.group())
        # print("match_z = ", match_z, match_z.group())

        NX = int(match_x.group(2))
        NY = int(match_y.group(2))
        NZ = int(match_z.group(2))
        box_size_found = True 
        print("Box size found, NX, NY, NZ = ", NX, NY, NZ)

    else: 
        print("All matches not found")

    return box_size_found, NX, NY, NZ


if __name__ == '__main__':

    basic = True  

    if basic:
    
        if len(sys.argv) >= 2:
            nprocs = int(sys.argv[1]) # number of threads to launch  
        else: 
            print("using default nprocs = 1")
            nprocs = 1 

        if len(sys.argv) >= 3:
            boundary_mesh_name = sys.argv[2]
        # compute masks for all 
        else:
            boundary_mesh_name = '2_aorta_remeshed_pt5mm_capped.vtp'

        # first make sure there is a times file 
        if not os.path.isfile('times.txt'):
            subprocess.call('visit -cli -nowin -s ~/heart_valves/scripts/write_times_file_visit.py', shell=True)

        times = []
        times_file = open('times.txt', 'r')
        for line in times_file:
            times.append(float(line)) 

        if times[0] == 0:
            dt = times[1]
        else: 
            dt = times[1] - times[0]

        # crop times for debug 
        # times = times[:5]

        point_data = True 
        
        box_size_found, NX, NY, NZ  = find_box_size()

        if not box_size_found:
            print("Log file not found. Using default box sizes.")

            NX=144
            NY=96
            NZ=224

            if "_192_" in os.getcwd():
                NX=96
                NY=48
                NZ=128

        print("NX, NY, NZ = ", NX, NY, NZ)

        basename = "eulerian_vars"
        nsteps = len(times)

        convert_to_point_data = True

        point_data = True 
        if point_data:
            label = '_restricted_points'
        else: 
            label = '_restricted_cells'

        extension = 'vtu'

        # grab this file if it's not here... 
        if not os.path.isfile(boundary_mesh_name):
            if os.path.isfile('~/heart_valves/' + boundary_mesh_name):
                shutil.copy('~/heart_valves/' + boundary_mesh_name, '.') 
            else: 
                raise FileNotFoundError("cannot find boundary_mesh_name file = ", boundary_mesh_name)


        boundary_mesh = pyvista.read(boundary_mesh_name)

        # selected = mesh.select_enclosed_points(boundary_mesh, tolerance=1.0e-10, inside_out=False, check_surface=True)
        # selected.clear_cell_arrays()

        run_all = True 

        if run_all:

            jobs = []
            for proc_num in range(nprocs):
                p = multiprocessing.Process(target=remove_eulerian_space, args=(basename, 
                                                                                boundary_mesh, 
                                                                                nsteps, 
                                                                                label, 
                                                                                extension,
                                                                                convert_to_point_data, 
                                                                                NX,
                                                                                NY, 
                                                                                NZ, 
                                                                                proc_num, 
                                                                                nprocs))
                jobs.append(p)
                p.start()

            for p in jobs:
                p.join()

            basename_out = basename + label

            write_pvd(basename_out, dt, nsteps, extension, nprocs_sim=1)



    debug_simple = False 
    if debug_simple:

        dir_name = 'eulerian_vars1365'    
        fname = 'eulerian_vars_combined_cells_1365.vtu'
        mesh = read_distributed_vtr(dir_name)    
        # mesh = mesh.cell_data_to_point_data()

        # mesh.save(fname)
        # mesh = pyvista.StructuredGrid(fname)

        fname_points = 'eulerian_vars_combined_point_1365.vtu'
        mesh_points = mesh.cell_data_to_point_data()
        mesh_points.save(fname_points)

        selected = mesh.select_enclosed_points(boundary_mesh, tolerance=1.0e-10, inside_out=False, check_surface=True)

        #breakpoint()

        mesh_inside = selected.threshold(0.5, scalars="SelectedPoints", all_scalars=False) 

        # mesh_inside = selected.threshold(0.5, all_scalars=True) 

        print("selected = ", selected)
        print("boundary_mesh = ", boundary_mesh)
        print("mesh_inside = ", mesh_inside)

        fname = 'eulerian_vars_cells_internal_1365.vtu'

        mesh_inside.save(fname)

        dargs = dict(show_edges=True)
        p = pyvista.Plotter()
        p.add_mesh(mesh_inside, color="Crimson", opacity=0.05, **dargs) 
        p.add_mesh(boundary_mesh, color="mintcream", opacity=1, **dargs) 
        p.show()


    aorta_single_frame = False
    if aorta_single_frame:

        frame_number = 421 

        point_data = True 
        NX=192
        NY=96
        NZ=256

        if "_192_" in os.getcwd():
            NX=96
            NY=48
            NZ=128

        basename = "eulerian_vars"

        if point_data:
            label = '_restricted_points'
        else: 
            label = '_restricted_cells'

        extension = 'vtu'

        # compute masks for all 
        boundary_mesh_name = '2_aorta_remeshed_pt5mm_capped.vtp'

        # grab this file if it's not here... 
        if not os.path.isfile(boundary_mesh_name):
            if os.path.isfile('~/heart_valves/' + boundary_mesh_name):
                shutil.copy('~/heart_valves/' + boundary_mesh_name, '.') 
            else: 
                raise FileNotFoundError("cannot find boundary_mesh_name file = ", boundary_mesh_name)


        boundary_mesh = pyvista.read(boundary_mesh_name)

        # selected = mesh.select_enclosed_points(boundary_mesh, tolerance=1.0e-10, inside_out=False, check_surface=True)
        # selected.clear_cell_arrays()
        remove_eulerian_space_single_frame(basename, 
                                            boundary_mesh, 
                                            frame_number, 
                                            label, 
                                            extension,
                                            point_data, 
                                            NX,
                                            NY, 
                                            NZ)




    pa_single_frame = False 
    if pa_single_frame:

        basename = "eulerian_vars"
        step = 3325
 
        label = '_restricted_cells'

        extension = 'vtu'

        dir_name = 'eulerian_vars0459'    
        fname = 'eulerian_vars_combined_cells_3325.vtu'
        mesh = read_distributed_vtr(dir_name)    

        mesh.save('eulerian_vars_combined_no_restrict_cells_3325.vtu')

        boundary_mesh_name = '3_1_three_end_crop_capped.stl'

        boundary_mesh = pyvista.read(boundary_mesh_name)

        selected = mesh.select_enclosed_points(boundary_mesh, tolerance=1.0e-10, inside_out=False, check_surface=True)

        mesh_inside = selected.threshold(0.5, scalars="SelectedPoints", all_scalars=True) 

        mesh_inside.save(fname)

        


