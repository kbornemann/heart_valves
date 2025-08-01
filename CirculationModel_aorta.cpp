// Filename: CirculationModel_aorta.cpp
// Created on 20 Aug 2007 by Boyce Griffith

// Modified 2019, Alexander D. Kaiser

#include "CirculationModel_aorta.h"
#include "pnpoly.h"
/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

#include <Eigen/Dense>
using namespace Eigen;

namespace
{
// Name of output file.
static const string DATA_FILE_NAME = "bc_data.m";

} 

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CirculationModel_aorta::CirculationModel_aorta(Pointer<Database> input_db, 
                                               const fourier_series_data *fourier_ventricle, 
                                               string ventricle_vertices_file_name,
                                               string aorta_vertices_file_name,
                                               const double  cycle_duration,
                                               const double  t_offset_bcs_unscaled, 
                                               const double  initial_time, 
                                               double P_initial_aorta,
                                               bool rcr_bcs_on,
                                               bool ventricle_0D_on, 
                                               bool P_initial_aorta_equal_to_ventricle,
                                               double rcr_on_time)
    : 
      d_object_name("circ_model_aorta"),  // constant name here  
      d_registered_for_restart(true),      // always true
      d_fourier_ventricle(fourier_ventricle), 
      d_cycle_duration(cycle_duration),
      d_t_offset_bcs_unscaled(t_offset_bcs_unscaled),
      d_current_idx_series(0),
      d_Q_ventricle(0.0), 
      d_Q_aorta(0.0),
      d_time(initial_time), 
      d_aorta_P_Wk(P_initial_aorta),
      d_p_extender_mean(0.0),
      d_p_extender_point(0.0),
      d_area_ventricle(0.0),
      d_area_aorta(0.0),
      d_area_initialized(false), 
      d_rcr_bcs_on(rcr_bcs_on), 
      d_ventricle_0D_on(ventricle_0D_on),
      d_P_initial_aorta_equal_to_ventricle(P_initial_aorta_equal_to_ventricle), 
      d_rcr_on_time(rcr_on_time)
{
    
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart databases.
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }
    
    if (d_rcr_bcs_on){
        if (input_db){
            // left and right equal for now
            d_aorta_R_proximal = input_db->getDouble("R_proximal");
            d_aorta_R_distal   = input_db->getDouble("R_distal");
            d_aorta_C          = input_db->getDouble("C");

            std::cout << "input db got values:\n";
            std::cout << "input db got values R_proximal = " << d_aorta_R_proximal << "\tR_distal = " << d_aorta_R_distal << "\tC = " << d_aorta_C << "\n";   
        }
        else {
            TBOX_ERROR("Must provide valid input_db");
        }
    }

    if (d_ventricle_0D_on){
        d_ventricle_0D = new ventricle_0D_model(input_db, d_cycle_duration, d_rcr_on_time); 
    } else{
        d_ventricle_0D = NULL; 
        if (d_fourier_ventricle == NULL){
            TBOX_ERROR("Must provide valid series for ventricle if d_ventricle_0D_on is off\n"); 
        }
    }

    // double x,x_prev,y,y_prev,z,z_prev; 
    double coord_normal, coord_normal_prev; 
    double tol = 1.0e-2; 

    // read vertices from file 
    ifstream ventricle_file(ventricle_vertices_file_name.c_str(), ios::in);

    if(!ventricle_file){
        TBOX_ERROR("Aorta file not found\n"); 
    }

    ventricle_file >> d_n_pts_ventricle; 
    

    double *ventricle_points_idx0 = new double[d_n_pts_ventricle];     
    double *ventricle_points_idx1 = new double[d_n_pts_ventricle]; 
    double *ventricle_points_idx2 = new double[d_n_pts_ventricle]; 

    int coords_flat[3] = {1, 1, 1};
    double x_prev, y_prev, z_prev;

    for (int i=0; i<d_n_pts_ventricle; i++){

        ventricle_file >> ventricle_points_idx0[i]; 
        ventricle_file >> ventricle_points_idx1[i]; 
        ventricle_file >> ventricle_points_idx2[i]; 

        if (i>0){
            if (fabs(x_prev - ventricle_points_idx0[i]) > tol){
                coords_flat[0] = 0;                
            }
            if (fabs(y_prev - ventricle_points_idx1[i]) > tol){
                coords_flat[1] = 0;                
            }
            if (fabs(z_prev - ventricle_points_idx2[i]) > tol){
                coords_flat[2] = 0;                
            }
        }

        x_prev = ventricle_points_idx0[i]; 
        y_prev = ventricle_points_idx1[i]; 
        z_prev = ventricle_points_idx2[i]; 

    }

    if (coords_flat[0]){
        
        if (coords_flat[1] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_ventricle_axis = 0;

        // checks whether in quarter of domain in normal direction 
        double tol_domain_x = (input_db->getDouble("X_HIGH") - input_db->getDouble("X_LOW"))/4.0;

        if (fabs(x_prev - input_db->getDouble("X_LOW")) < tol_domain_x){
            d_ventricle_side = 0;
        }
        else if (fabs(x_prev - input_db->getDouble("X_HIGH")) < tol_domain_x){
            d_ventricle_side = 1;
        }
        else{
            TBOX_ERROR("Flat coodidinate not near boundary");
        }

        // gets y and z coords         
        delete[] ventricle_points_idx0;
        d_ventricle_points_idx1 = ventricle_points_idx1;
        d_ventricle_points_idx2 = ventricle_points_idx2;

    }
    else if (coords_flat[1]){
        if (coords_flat[0] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_ventricle_axis = 1;        

        // checks whether in quarter of domain in normal direction 
        double tol_domain_y = (input_db->getDouble("Y_HIGH") - input_db->getDouble("Y_LOW"))/4.0;

        if (fabs(y_prev - input_db->getDouble("Y_LOW")) < tol_domain_y){
            d_ventricle_side = 0;
        }
        else if (fabs(y_prev - input_db->getDouble("Y_HIGH")) < tol_domain_y){
            d_ventricle_side = 1;
        }
        else{
            TBOX_ERROR("Flat coodidinate not near boundary");
        }

        // gets x and z coords         
        d_ventricle_points_idx1 = ventricle_points_idx0;
        delete[] ventricle_points_idx1;
        d_ventricle_points_idx2 = ventricle_points_idx2;        

    }
    else if (coords_flat[2]){
        if (coords_flat[0] || coords_flat[1]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_ventricle_axis = 2;        

        // checks whether in quarter of domain in normal direction 
        double tol_domain_z = (input_db->getDouble("Z_HIGH") - input_db->getDouble("Z_LOW"))/4.0;
        
        if (fabs(z_prev - input_db->getDouble("Z_LOW")) < tol_domain_z){
            d_ventricle_side = 0;
        }
        else if (fabs(z_prev - input_db->getDouble("Z_HIGH")) < tol_domain_z){
            d_ventricle_side = 1;
        }
        else{
            TBOX_ERROR("Flat coodidinate not near boundary");
        }

        // gets x and y coords         
        d_ventricle_points_idx1 = ventricle_points_idx0;
        d_ventricle_points_idx2 = ventricle_points_idx1;
        delete[] ventricle_points_idx2;

    }
    else{
        TBOX_ERROR("Flat coodidinate not found in boundary condition");
    }

    pout << "to ventricle file close\n"; 
    pout << "found d_ventricle_side = " << d_ventricle_side << ", d_ventricle_axis = " << d_ventricle_axis << "\n";
    ventricle_file.close(); 



    // read vertices from file 
    ifstream aorta_file(aorta_vertices_file_name.c_str(), ios::in);

    if(!aorta_file){
        TBOX_ERROR("Aorta file not found\n"); 
    }

    aorta_file >> d_n_pts_aorta; 
    
    double *aorta_points_idx0 = new double[d_n_pts_aorta];     
    double *aorta_points_idx1 = new double[d_n_pts_aorta]; 
    double *aorta_points_idx2 = new double[d_n_pts_aorta]; 

    coords_flat[0] = 1; 
    coords_flat[1] = 1; 
    coords_flat[2] = 1; 

    for (int i=0; i<d_n_pts_aorta; i++){

        aorta_file >> aorta_points_idx0[i]; 
        aorta_file >> aorta_points_idx1[i]; 
        aorta_file >> aorta_points_idx2[i]; 

        if (i>0){
            if (fabs(x_prev - aorta_points_idx0[i]) > tol){
                coords_flat[0] = 0;                
            }
            if (fabs(y_prev - aorta_points_idx1[i]) > tol){
                coords_flat[1] = 0;                
            }
            if (fabs(z_prev - aorta_points_idx2[i]) > tol){
                coords_flat[2] = 0;                
            }
        }

        x_prev = aorta_points_idx0[i]; 
        y_prev = aorta_points_idx1[i]; 
        z_prev = aorta_points_idx2[i]; 

    }

    if (coords_flat[0]){
        
        if (coords_flat[1] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_aorta_axis = 0;

        if (fabs(x_prev - input_db->getDouble("X_LOW")) < tol){
            d_aorta_side = 0;
        }
        else if (fabs(x_prev - input_db->getDouble("X_HIGH")) < tol){
            d_aorta_side = 1;
        }
        else{
            TBOX_ERROR("Flat coodidinate not near boundary");
        }

        // gets y and z coords         
        delete[] aorta_points_idx0;
        d_aorta_points_idx1 = aorta_points_idx1;
        d_aorta_points_idx2 = aorta_points_idx2;

    }
    else if (coords_flat[1]){
        if (coords_flat[0] || coords_flat[2]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_aorta_axis = 1;        

        if (fabs(y_prev - input_db->getDouble("Y_LOW")) < tol){
            d_aorta_side = 0;
        }
        else if (fabs(y_prev - input_db->getDouble("Y_HIGH")) < tol){
            d_aorta_side = 1;
        }
        else{
            TBOX_ERROR("Flat coodidinate not near boundary");
        }

        // gets x and z coords         
        d_aorta_points_idx1 = aorta_points_idx0;
        delete[] aorta_points_idx1;
        d_aorta_points_idx2 = aorta_points_idx2;        

    }
    else if (coords_flat[2]){
        if (coords_flat[0] || coords_flat[1]){
            TBOX_ERROR("More than one coordinate is flat");
        }

        d_aorta_axis = 2;        
        
        if (fabs(z_prev - input_db->getDouble("Z_LOW")) < tol){
            d_aorta_side = 0;
        }
        else if (fabs(z_prev - input_db->getDouble("Z_HIGH")) < tol){
            d_aorta_side = 1;
        }
        else{
            TBOX_ERROR("Flat coodidinate not near boundary");
        }

        // gets x and y coords         
        d_aorta_points_idx1 = aorta_points_idx0;
        d_aorta_points_idx2 = aorta_points_idx1;
        delete[] aorta_points_idx2;

    }
    else{
        TBOX_ERROR("Flat coodidinate not found in boundary condition");
    }

    pout << "to aorta file close\n"; 
    pout << "found d_aorta_side = " << d_aorta_side << ", d_aorta_axis = " << d_aorta_axis << "\n";

    pout << "to aorta file close\n"; 
    aorta_file.close(); 


    // misc scalars 
    if (!from_restart){

        // initial ventricle pressure 
        if (d_ventricle_0D_on){
            d_ventricle_P = d_ventricle_0D->d_P_ventricle;
        }
        else {
            d_ventricle_P = MMHG_TO_CGS * d_fourier_ventricle->values[0];
        }

        // initial aorta pressure 
        if (d_P_initial_aorta_equal_to_ventricle){
            d_aorta_P = d_ventricle_P; 
        }
        else{
            d_aorta_P = P_initial_aorta; 
        }

        // temp value unused 
        d_P_min_linear_interp = 0.0; 
    }

    d_p_equal_fraction = 0.1; 

    /* 
    // get fourier series value for pressure interpolation 
    
    double t_reduced = d_p_equal_fraction * d_rcr_on_time;
    double t_scaled = t_reduced * (d_fourier_ventricle->L  / d_cycle_duration);
    double t_scaled_offset = t_scaled + d_t_offset_bcs_unscaled;
    unsigned int k = (unsigned int) floor(t_scaled_offset / (d_fourier_ventricle->dt));
    unsigned int idx = k % (d_fourier_ventricle->N_times);

    // pressure value half way through initialization period 
    d_P_min_linear_interp = MMHG_TO_CGS * d_fourier_ventricle->values[idx]; 
    */ 

    pout << "passed contstructor\n"; 

    pout << "initial aorta pressure = " << P_initial_aorta << ", P_wk = " << d_aorta_P << "\n"; 

    return;
} // CirculationModel

CirculationModel_aorta::~CirculationModel_aorta()
{
    return;
} // ~CirculationModel_aorta


void CirculationModel_aorta::advanceTimeDependentData(const double dt,
                                                        const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                        const int U_idx,
                                                        const int /*P_idx*/,
                                                        const int /*wgt_cc_idx*/,
                                                        const int wgt_sc_idx)
{
    // Compute the mean flow rates in the vicinity of the inflow and outflow
    // boundaries.
    
    double Q_ventricle_local = 0.0; 
    double Q_aorta_local = 0.0; 

    double area_ventricle_local = 0.0; 
    double area_aorta_local = 0.0; 

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            if (pgeom->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(U_idx);
                Pointer<SideData<NDIM, double> > wgt_sc_data = patch->getPatchData(wgt_sc_idx);
                const Box<NDIM>& patch_box = patch->getBox();
                const double* const x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                double dV = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    dV *= dx[d];
                }

                for(int axis=0; axis<3; axis++)
                {
                    for(int side=0; side<2; side++)
                    {
                        const bool is_lower = (side == 0);
                        if (pgeom->getTouchesRegularBoundary(axis, side))
                        {
                            
                            IBTK::Vector n;
                            for (int d = 0; d < NDIM; ++d)
                            {
                                n[d] = axis == d ? (is_lower ? -1.0 : +1.0) : 0.0;
                            }
                            Box<NDIM> side_box = patch_box;
                            if (is_lower)
                            {
                                side_box.lower(axis) = patch_box.lower(axis);
                                side_box.upper(axis) = patch_box.lower(axis);
                            }
                            else
                            {
                                side_box.lower(axis) = patch_box.upper(axis) + 1;
                                side_box.upper(axis) = patch_box.upper(axis) + 1;
                            }
                            for (Box<NDIM>::Iterator b(side_box); b; b++)
                            {
                                const SAMRAI::hier::Index<NDIM>& i = b();

                                double X[NDIM];
                                for (int d = 0; d < NDIM; ++d)
                                {
                                    X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == axis ? 0.0 : 0.5));
                                }

                                double X_in_plane_1 = 0.0; 
                                double X_in_plane_2 = 0.0; 
                                if (axis == 0)
                                {
                                    X_in_plane_1 = X[1]; 
                                    X_in_plane_2 = X[2]; 
                                }
                                else if (axis == 1)
                                {
                                    X_in_plane_1 = X[0]; 
                                    X_in_plane_2 = X[2]; 
                                }
                                else if (axis == 2)
                                {
                                    X_in_plane_1 = X[0]; 
                                    X_in_plane_2 = X[1]; 
                                }
                                else{
                                    TBOX_ERROR("Invalid value of axis\n"); 
                                }

                                const int in_ventricle  = this->point_in_ventricle(X_in_plane_1, X_in_plane_2, axis, side);
                                const int in_aorta      = this->point_in_aorta    (X_in_plane_1, X_in_plane_2, axis, side);

                                if (in_ventricle && in_aorta){
                                    TBOX_ERROR("Position is within two inlets and outlets, should be impossible\n"); 
                                }

                                if (in_ventricle)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_ventricle_local += (*U_data)(i_s)* n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_ventricle_local += dA;
                                        }

                                    }
                                }

                                if (in_aorta)
                                {
                                    const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                    if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                    {
                                        double dA = dV / dx[axis];
                                        Q_aorta_local += (*U_data)(i_s) * n[axis] * dA;

                                        if (!d_area_initialized){
                                            area_aorta_local += dA;
                                        }

                                    }
                                }

                            }
                        }
                    }
                }
            }
        }
    }

    d_Q_ventricle = SAMRAI_MPI::sumReduction(Q_ventricle_local);
    d_Q_aorta     = SAMRAI_MPI::sumReduction(Q_aorta_local);

    if (!d_area_initialized){
        d_area_ventricle   = SAMRAI_MPI::sumReduction(area_ventricle_local);
        d_area_aorta       = SAMRAI_MPI::sumReduction(area_aorta_local);  
        d_area_initialized = true;       
    }

    if (d_rcr_bcs_on){
        
        // The downstream pressure is determined by a three-element Windkessel model.
        
        if ((d_P_initial_aorta_equal_to_ventricle) && (d_time < (d_p_equal_fraction*d_rcr_on_time))){
            // equal to ventricle for half the time 
            d_aorta_P = d_ventricle_P; 
            d_P_min_linear_interp = d_ventricle_P; 
        }
        else if ((d_P_initial_aorta_equal_to_ventricle) && (d_time < d_rcr_on_time)){
            // linear interpolation to pressurize 

            // linear interpolation to pressurize 
            // resistance hooked to linearly interpolated pressures 
            // distal pressure is linearly interpolated 
            double P_distal_temp = ((d_time  -     d_rcr_on_time) / (d_p_equal_fraction * d_rcr_on_time - d_rcr_on_time)) * d_P_min_linear_interp + 
                                   ((d_time  - d_p_equal_fraction * d_rcr_on_time) / (d_rcr_on_time - d_p_equal_fraction * d_rcr_on_time)) * d_aorta_P_Wk; // wk pressure is the end pressure for the interpolation

            // then gets the proximal resistance, no capacitor 
            d_aorta_P = P_distal_temp + d_aorta_R_proximal * d_Q_aorta; 

        }
        else{
            d_aorta_P_Wk = ((d_aorta_C / dt) * d_aorta_P_Wk + d_Q_aorta) / (d_aorta_C / dt + 1.0 / d_aorta_R_distal);        
            d_aorta_P = d_aorta_P_Wk + d_aorta_R_proximal * d_Q_aorta;
        }
    }

    // print_summary();

    // bool debug_out_areas = false; 
    // if (debug_out_areas){
    //     pout << "d_area_ventricle = " << d_area_ventricle << "\n"; 
    //     pout << "d_area_aorta = " << d_area_aorta << "\n"; 
    // }

    d_time += dt; 

    if (d_ventricle_0D_on){

        // d_Q_ventricle is outflow according to surface normal 
        // actual Q_ventricle gets negative
        d_ventricle_0D->advanceTimeDependentData(dt, d_time, -d_Q_ventricle);
        // local ventricle pressure copied from lv circ model 
        d_ventricle_P = d_ventricle_0D->d_P_ventricle; 

    }
    else{
        // compute which index in the Fourier series we need here 
        // always use a time in current cycle 
        double t_reduced = d_time - d_cycle_duration * floor(d_time/d_cycle_duration); 

        // fourier series has its own period, scale to that 
        double t_scaled = t_reduced * (d_fourier_ventricle->L  / d_cycle_duration); 

        // start offset some arbitrary time in the cardiac cycle, but this is relative to the series length 
        double t_scaled_offset = t_scaled + d_t_offset_bcs_unscaled; 

        // Fourier data here
        // index without periodicity 
        unsigned int k = (unsigned int) floor(t_scaled_offset / (d_fourier_ventricle->dt));
        
        // // take periodic reduction
        d_current_idx_series = k % (d_fourier_ventricle->N_times);

        d_ventricle_P = MMHG_TO_CGS * d_fourier_ventricle->values[d_current_idx_series];

        // bool debug_out = false; 
        // if (debug_out){
        //     pout << "circ mode: d_time = " << d_time << ", d_current_idx_series = " << d_current_idx_series << "\n"; 
        //     pout << "t_reduced = " << t_reduced << " t_scaled = " << t_scaled << " t_scaled_offset = " << t_scaled_offset << "\n"; 
        //     pout << "k (unreduced idx) = " << k << " d_current_idx_series = " << d_current_idx_series << "\n\n"; 
        // }        
    }




    writeDataFile(); 

} // advanceTimeDependentData

void CirculationModel_aorta::set_Q_valve(double Q_valve){
    d_Q_valve = Q_valve; 
}

void CirculationModel_aorta::set_extender_pressures(double p_extender_mean, double p_extender_point){
    d_p_extender_mean  = p_extender_mean; 
    d_p_extender_point = p_extender_point;
}



void
CirculationModel_aorta::putToDatabase(Pointer<Database> db)
{

    db->putInteger("d_current_idx_series", d_current_idx_series); 
    db->putDouble("d_ventricle_P", d_ventricle_P); 
    db->putDouble("d_Q_ventricle", d_Q_ventricle); 
    db->putDouble("d_Q_aorta", d_Q_aorta);
    db->putDouble("d_Q_valve", d_Q_valve);
    db->putDouble("d_aorta_P", d_aorta_P);
    db->putDouble("d_aorta_P_Wk", d_aorta_P_Wk);
    db->putDouble("d_p_extender_mean", d_p_extender_mean);
    db->putDouble("d_p_extender_point", d_p_extender_point);
    db->putDouble("d_time", d_time); 
    db->putBool("d_rcr_bcs_on", d_rcr_bcs_on); 
    db->putBool("d_ventricle_0D_on", d_ventricle_0D_on); 
    db->putDouble("d_P_min_linear_interp", d_P_min_linear_interp);
    return; 
} // putToDatabase

void CirculationModel_aorta::print_summary(){

    double P_ventricle = d_ventricle_P / MMHG_TO_CGS; 
    double P_aorta = d_aorta_P / MMHG_TO_CGS;

    if (!d_rcr_bcs_on){
        TBOX_ERROR("Not implemented\n"); 
        // P_aorta        = d_fourier_aorta->values[d_current_idx_series];
    }

    pout << "rcr_bcs_on = " << d_rcr_bcs_on << "\n"; 
    pout << "% time \t       P_ventricle (mmHg)\t   P_aorta (mmHg)\t  Q_ventricle (ml/s)\t    d_Q_aorta (ml/s)\t  d_Q_valve (ml/s)\t  Q_current_idx_series \t idx" ;
    pout << "\t aorta_P_Wk \t p_extender_mean \t p_extender_point "; 
    pout << "\n";
    pout << d_time << " " << P_ventricle <<  " " << P_aorta << " " << d_Q_ventricle << " " << d_Q_aorta << " " << d_Q_valve << " " << d_current_idx_series; 
    pout  << " " << d_aorta_P_Wk << " " << d_p_extender_mean/MMHG_TO_CGS << " " << d_p_extender_point/MMHG_TO_CGS; 
    pout << "\n";

}

void CirculationModel_aorta::print_bc_debug(){

    pout << "d_ventricle_side = " << d_ventricle_side << ", d_ventricle_axis = " << d_ventricle_axis << "\n"; 

    pout << "ventricle_points:\n";
    for(int i=0; i<d_n_pts_ventricle; i++){
        pout << d_ventricle_points_idx1[i] << ", " << d_ventricle_points_idx2[i] << "\n";
    }
    pout << "\n";

    pout << "d_aorta_side = " << d_aorta_side << ", d_aorta_axis = " << d_aorta_axis << "\n"; 
    pout << "aorta_points:\n";
    for(int i=0; i<d_n_pts_aorta; i++){
        pout << d_aorta_points_idx1[i] << ", " << d_aorta_points_idx2[i] << "\n";
    }
    pout << "\n";    

}

int CirculationModel_aorta::point_in_ventricle(double testx, double testy, int axis, int side){
    // checks whether given point is in right ventricle

    // quick exit for correct side and axis 
    if ((axis != d_ventricle_axis) || (side != d_ventricle_side))
        return 0; 

    return pnpoly(d_n_pts_ventricle, d_ventricle_points_idx1, d_ventricle_points_idx2, testx, testy); 
}

int CirculationModel_aorta::point_in_aorta(double testx, double testy, int axis, int side){
    // checks whether given point is in right ventricle

    // quick exit for correct side and axis 
    if ((axis != d_aorta_axis) || (side != d_aorta_side))
        return 0; 

    return pnpoly(d_n_pts_aorta, d_aorta_points_idx1, d_aorta_points_idx2, testx, testy); 
}


void CirculationModel_aorta::write_plot_code()
{

    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);
        fout.setf(ios_base::scientific);
        fout.setf(ios_base::showpos);
        fout.precision(10);
        fout << "];\n";  
        fout << "fig = figure;\n";  
        fout << "times            =  bc_vals(:,1);\n"; 
        fout << "p_lv             =  bc_vals(:,2);\n"; 
        fout << "p_aorta          =  bc_vals(:,3); \n"; 
        fout << "q_ventricle      = -bc_vals(:,4);\n"; 
        fout << "q_aorta          =  bc_vals(:,5);\n"; 
        fout << "q_valve          =  bc_vals(:,6);\n"; 
        fout << "p_wk             =  bc_vals(:,7);\n"; 
        fout << "p_extender_mean  =  bc_vals(:,8);\n"; 
        fout << "p_extender_point =  bc_vals(:,9);\n"; 
        fout << "subplot(2,1,1)\n"; 
        fout << "plot(times, p_aorta, 'k')\n"; 
        fout << "hold on\n"; 
        fout << "plot(times, p_wk, ':k')\n"; 
        fout << "plot(times, p_lv, '--k')\n"; 
        fout << "%plot(times, p_extender_mean)\n"; 
        fout << "plot(times, p_extender_point)\n";                         
        fout << "%legend('P_{Ao}', 'P_{Wk}', 'P_{LV}', 'P extender mean', 'P extender point', 'Location','NorthEastOutside');\n"; 
        fout << "legend('P_{Ao}', 'P_{Wk}', 'P_{LV}', 'P extender point', 'Location','NorthEastOutside');\n"; 
        fout << "xlabel('t (s)');\n"; 
        fout << "ylabel('P (mmHg)');\n"; 
        fout << "subplot(2,1,2)\n"; 
        fout << "plot(times, q_aorta, 'k')\n"; 
        fout << "hold on\n"; 
        fout << "dt = times(2,1) - times(1);\n"; 
        fout << "net_flux = dt*cumsum(q_aorta);\n"; 
        fout << "plot(times, net_flux, '--k')\n"; 
        fout << "plot(times, q_ventricle)\n"; 
        fout << "% plot(times, q_valve)\n"; 
        fout << "plot(bc_vals(:,1), 0*net_flux, ':k')\n"; 
        fout << "%legend('Q', 'net Q', 'Q ventricle', 'Q valve', 'Location','NorthEastOutside')\n"; 
        fout << "legend('Q', 'net Q', 'Q ventricle', 'Location','NorthEastOutside')\n"; 
        fout << "xlabel('t (s)')\n"; 
        fout << "ylabel('Flow (ml/s), Net Flow (ml)')\n"; 
        fout << "set(fig, 'Position', [100, 100, 1000, 750])\n"; 
        fout << "set(fig,'PaperPositionMode','auto')\n"; 
        fout << "printfig(fig, 'bc_model_variables')\n"; 
        fout << "diary notes.txt\n"; 
        fout << "min_p_aorta_after_first_beat = min(p_aorta(floor(end/3):end))\n"; 
        fout << "max_p_aorta_after_first_beat = max(p_aorta(floor(end/3):end))\n"; 
        fout << "mean_p_aorta = mean(p_aorta)\n"; 
        fout << "mean_p_wk    = mean(p_wk)\n"; 
        fout << "max_p_exender_mean = max(p_extender_mean)\n"; 
        fout << "max_p_exender_point = max(p_extender_point)\n"; 
        fout << "mean_p_lv    = mean(p_lv)\n"; 
        fout << "end_flux    = net_flux(end)\n"; 
        fout << "max_flux    = max(net_flux)\n"; 
        fout << "diary off\n"; 
        fout << "save bc_data.mat\n";

    }
    return;

}



/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
    CirculationModel_aorta::writeDataFile() const
{
    static const int mpi_root = 0;
    if (SAMRAI_MPI::getRank() == mpi_root)
    {
        static bool file_initialized = false;
        const bool from_restart = RestartManager::getManager()->isFromRestart();
        if (!from_restart && !file_initialized)
        {
            ofstream fout(DATA_FILE_NAME.c_str(), ios::out);
            fout << "% time \t P_ventricle (mmHg)\t P_aorta (mmHg)\t d_Q_ventricle (ml/s)\t d_Q_aorta (ml/s) \t d_Q_valve (ml/s) \t " 
                 << "d_aorta_P_Wk (mmHg) \t d_p_extender_mean (mmHg) \t d_p_extender_point (mmHg)"; 

            if (d_ventricle_0D_on){
                fout << "\t d_V_ventricle \t d_V_rest_ventricle \t d_Elas \t d_act_temp \t d_Q_in \t d_P_ventricle_in \t d_P_lvot_upstream"; 
            }

            fout << "\n" 
                 << "bc_vals = [";
            file_initialized = true;
        }

        ofstream fout(DATA_FILE_NAME.c_str(), ios::app);

        fout.setf(ios_base::scientific);
        fout.setf(ios_base::showpos);
        fout.precision(10);

        fout << d_time;

        double P_ventricle = d_ventricle_P / MMHG_TO_CGS; 

        double P_aorta = 0.0; 

        if (d_rcr_bcs_on){
            P_aorta        = d_aorta_P/MMHG_TO_CGS;
        }
        else{
            TBOX_ERROR("not implemented\n"); 
            // P_aorta        = d_fourier_aorta->values[d_current_idx_series];
        }

        fout << " " << P_ventricle <<  " " << P_aorta;
        fout << " " << d_Q_ventricle << " " << d_Q_aorta << " " << d_Q_valve;         
        fout << " " << d_aorta_P_Wk/MMHG_TO_CGS;
        fout << " " << d_p_extender_mean/MMHG_TO_CGS;
        fout << " " << d_p_extender_point/MMHG_TO_CGS;    

        if (d_ventricle_0D_on){
            fout << " " << d_ventricle_0D->d_V_ventricle;
            fout << " " << d_ventricle_0D->d_V_rest_ventricle; 
            fout << " " << d_ventricle_0D->d_Elas;
            fout << " " << d_ventricle_0D->d_act_temp;
            fout << " " << d_ventricle_0D->d_Q_in;            
            fout << " " << d_ventricle_0D->d_P_ventricle_in/MMHG_TO_CGS;
            fout << " " << d_ventricle_0D->d_P_lvot_upstream/MMHG_TO_CGS;
        }

        fout << "; \n";

    }

    return;
} // writeDataFile

void
CirculationModel_aorta::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to " << d_object_name << " not found in restart file.");
    }

    d_current_idx_series  = db->getInteger("d_current_idx_series"); 
    d_ventricle_P         = db->getDouble("d_ventricle_P");
    d_Q_ventricle         = db->getDouble("d_Q_ventricle"); 
    d_Q_aorta             = db->getDouble("d_Q_aorta");
    d_Q_valve             = db->getDouble("d_Q_valve");
    d_aorta_P             = db->getDouble("d_aorta_P");
    d_aorta_P_Wk          = db->getDouble("d_aorta_P_Wk");
    d_p_extender_mean     = db->getDouble("d_p_extender_mean");
    d_p_extender_point    = db->getDouble("d_p_extender_point");
    d_time                = db->getDouble("d_time");
    d_rcr_bcs_on          = db->getBool("d_rcr_bcs_on"); 
    d_ventricle_0D_on     = db->getBool("d_ventricle_0D_on"); 
    d_P_min_linear_interp = db->getDouble("d_P_min_linear_interp");

    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////