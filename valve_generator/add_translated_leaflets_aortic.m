function valve = add_translated_leaflets_aortic(valve)

N_leaflets = valve.leaflets(1).N_leaflets; 

if length(valve.leaflets) ~= 1 
    error('rotation implmented only with one initial leaflet'); 
end 

if N_leaflets == 1
    warning('No rotation of one single leaflet')
    return; 
end 

th = 2*pi/N_leaflets; 

k_max = valve.leaflets(1).k_max; 

for n = 2:N_leaflets
    
    valve.leaflets(n) = valve.leaflets(1);
    
    % rotate all cols 
    for k = 1:k_max
        valve.leaflets(n).X(:,:,k) = rotation_matrix_z((n-1)*th) * valve.leaflets(n).X(:,:,k);
    end 
    
end

f = fopen('aortic_annulus_truncal_postop.vertex', 'r'); 
vertices_ring_bdry = fscanf(f, '%f'); 
fclose(f); 

% first is number of vertices 
n_pts_ring_from_file = vertices_ring_bdry(1); 

% crop to keep just this 
vertices_ring_bdry = vertices_ring_bdry(2:end); 

vertices_ring_bdry = reshape(vertices_ring_bdry, 3, []); 
[~, n_pts_ring] = size(vertices_ring_bdry);

