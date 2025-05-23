function fig = dissection_plot_rest_height_aortic(valve, fig)

if ~exist('fig', 'var')
    fig = figure; 
end 

hold off 

if ~isfield(valve, 'name') || ~strcmp(valve.name, 'aortic') 
    error('dissection_plot_rest_height_aortic called with non aortic type'); 
    return; 
end 

% loaded (current configuration) lengths
leaflet = valve.leaflets(1); 
X_current              = leaflet.X; 
j_max                  = leaflet.j_max; 
k_max                  = leaflet.k_max; 
is_internal            = leaflet.is_internal; 
is_bc                  = leaflet.is_bc; 
N_each                 = leaflet.N_each; 

center_leaflet_idx = N_each/2 + 1; 
j = center_leaflet_idx; 
radial_leaflet_height_loaded = 0; 
for k=1:(k_max-1)
    if is_internal(j,k) || is_bc(j,k)
        j_nbr_tmp = j; 
        k_nbr_tmp = k + 1; 

        [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
        if ~valid 
            error('trying to compute lengths with an invalid rest length')
        end
        
        X = X_current(:,j,k);
        X_nbr = X_current(:,j_nbr,k_nbr); 
        
        radial_leaflet_height_loaded = radial_leaflet_height_loaded + norm(X - X_nbr);
        
    end 
end 

free_edge_length_single_loaded = 0; 
for j=1:N_each
    k=k_max; 
    
    j_nbr_tmp = j+1; 
    k_nbr_tmp = k; 
    [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
    if ~valid 
        error('trying to compute lengths with an invalid rest length')
    end
    
    X = X_current(:,j,k);
    X_nbr = X_current(:,j_nbr,k_nbr); 
    
    free_edge_length_single_loaded = free_edge_length_single_loaded + norm(X - X_nbr);
end 


valve_with_reference = valve; 
valve_with_reference = rmfield(valve_with_reference, 'leaflets'); 
valve_with_reference.leaflets(1) = set_rest_lengths_and_constants_aortic(valve.leaflets(1), valve); 

leaflet = valve_with_reference.leaflets(1); 

r_commissure = valve.skeleton.r_commissure; 
radius = valve.skeleton.r; 

X_current              = leaflet.X; 
j_max                  = leaflet.j_max; 
k_max                  = leaflet.k_max; 
is_internal            = leaflet.is_internal; 
is_bc                  = leaflet.is_bc; 
R_v                    = leaflet.R_v; 
R_u                    = leaflet.R_u; 
N_each                 = leaflet.N_each; 
N_leaflets             = leaflet.N_leaflets; 

circ_cumulative_sum = linspace(0, 2*pi*radius / N_leaflets, j_max); 

positions_x = nan(j_max,k_max); 
positions_y = nan(j_max,k_max); 

center_leaflet_idx = N_each/2 + 1; 

for j=1:j_max
    for k=1:k_max
        if is_internal(j,k) || is_bc(j,k)
            positions_x(j,k) = circ_cumulative_sum(j); 
            
            positions_y(j,k) = X_current(3,j,1); % z component at bottom 

            if j == center_leaflet_idx 
                'stop'; 
            end 
            
            if k > 1 
                for k_tmp = 1:(k-1)
                    
                    j_nbr_tmp = j; 
                    k_nbr_tmp = k_tmp + 1; 
                    
                    [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k_tmp, j_nbr_tmp, k_nbr_tmp); 

                    if ~valid 
                        error('trying to compute lengths with an invalid rest length')
                    end 
                    positions_y(j,k) = positions_y(j,k) + R_v(j_spr,k_spr); 

                end 
                
            end 
            
        end 
    end 
end 

subplot(2,1,1)
plot(positions_x, positions_y, 'k.'); 
axis equal
xlabel('circumferential (cm)')
ylabel('radial (cm)')
title('aortic, radial height only')



length_one_leaflet_free_edge = zeros(k_max,1); 
for j=1:N_each
    for k=1:k_max
        if is_internal(j,k) || is_bc(j,k)
            j_nbr_tmp = j+1; 
            k_nbr_tmp = k; 
            [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 

            if ~valid 
                error('trying to compute lengths with an invalid rest length')
            end

            length_one_leaflet_free_edge(k) = length_one_leaflet_free_edge(k) + R_u(j_spr,k_spr); 
        end 
    end 
end 

subplot(2,1,2)
j = center_leaflet_idx; 
for k=1:k_max
    y = [positions_y(j,k) positions_y(j,k)]; 
    x = length_one_leaflet_free_edge(k)/2 * [-1, 1]; 
    plot(x,y,'k'); 
    hold on 
end 
    
radial_leaflet_height = positions_y(center_leaflet_idx,k_max) - positions_y(center_leaflet_idx,1)

% figure; 
height_percentage_ref_aortic = cumsum(R_v(center_leaflet_idx,:)) / sum(R_v(center_leaflet_idx,:)); 
% plot(1:k_max,height_percentage_ref_aortic)
% title('percent of leaflet height')
% xlabel('k')
% ylabel('height so far (fraction)')
% 
% figure; 
height_value_ref_aortic = cumsum(R_v(center_leaflet_idx,:)); 
% plot(1:k_max,height_value_ref_aortic)
% title('absolute leaflet height in reference config')
% xlabel('k')
% ylabel('height so far (absolute)')

coaptation_height_goal = 0.17*2*radius 
coaptation_zone_start_index = find(height_value_ref_aortic > (radial_leaflet_height - coaptation_height_goal),1)
coaptation_zone_start_fraction = coaptation_zone_start_index/k_max 

% r = d/2 (clearly)
% h = 0.71d = 1.42r, height to commissure 
% c = .17d = .34r, coaptation height
% ?f = .62 = 1.24r, free edge(is this entire free edge length?). this must be half length. 
% lc = .7d = 1.4r. the inner length of the dotted line on the left,  leaflet height in 2d. 

stretch_circ = valve.strain_circ + 1; 
stretch_rad = valve.strain_rad + 1; 

if isfield(leaflet, 'variety') && strcmp(leaflet.variety, 'bicuspid') 
    fprintf('Rest length height summary, true bicuspid:\n'); 
    fprintf("Radius                             = %f\n", radius); 
    fprintf("One half vbr circumference         = %f\t", pi*radius); 
    fprintf("One half inter comm circumference  = %f\n", pi*r_commissure); 
    fprintf("Geometric height goal (1.3 r)      = %f\n", 1.3*radius);
    fprintf("Circ free edge loaded length       = %f\t", free_edge_length_single_loaded); 
    fprintf("Circ free edge rest length         = %f\n", length_one_leaflet_free_edge(k_max));
    fprintf("Ratio loaded free edge to vbr d    = %f\t", free_edge_length_single_loaded/(2*radius));     
    fprintf("Ratio rest   free edge to vbr d    = %f\n", length_one_leaflet_free_edge(k_max)/(2*radius));   
    fprintf("Ratio loaded free edge to stj d    = %f\t", free_edge_length_single_loaded/(2*r_commissure));     
    fprintf("Ratio rest   free edge to stj d    = %f\n", length_one_leaflet_free_edge(k_max)/(2*r_commissure));    
    fprintf("Radial height loaded length        = %f\t", radial_leaflet_height_loaded); 
    fprintf("Radial height rest length          = %f\n", radial_leaflet_height); 
    fprintf("Ratio loaded radial h to vbr d     = %f\t", radial_leaflet_height_loaded/(2*radius));     
    fprintf("Ratio rest   radial h to vbr d     = %f\n", radial_leaflet_height/(2*radius));   
    fprintf("Ratio loaded radial h to stj d     = %f\t", radial_leaflet_height_loaded/(2*r_commissure));     
    fprintf("Ratio rest   radial h to stj d     = %f\n", radial_leaflet_height/(2*r_commissure));    
    
    
else 
    fprintf('Rest length height summary:\n'); 
    fprintf("And targets based on Swanson and Clark 1973\n")
    fprintf("Radius                       = %f\n", radius); 
    fprintf("Coaptation height            = %f\n", 0.17*2*radius); 
    fprintf("One third circumference      = %f\n", 2*pi*radius/3); 
    fprintf("Free edge target             = %f\t", 2*1.24*radius); 
    fprintf("Free edge rest target        = %f\n", 2*1.24*radius / stretch_circ); 
    fprintf("Circ free edge loaded length = %f\t", free_edge_length_single_loaded); 
    fprintf("Circ free edge rest length   = %f\n", length_one_leaflet_free_edge(k_max));
    fprintf("Leaflet height target        = %f\t", 1.4*radius); 
    fprintf("Leaflet height rest target   = %f\n", 1.4*radius / stretch_rad); 
    fprintf("Radial height loaded length  = %f\t", radial_leaflet_height_loaded); 
    fprintf("Radial height rest length    = %f\n", radial_leaflet_height); 
end 

% 
% j_center_anterior  = leaflet.j_range_anterior(floor(end/2)); 
% j_center_posterior = leaflet.j_range_posterior(floor(end/2)); 
% 
% k_center_anterior = find(~isnan(positions_y(j_center_anterior,:)), 1); 
% height_anterior = abs(positions_y(j_center_anterior, k_center_anterior)); 
% 
% k_center_posterior = find(~isnan(positions_y(j_center_posterior,:)), 1); 
% height_posterior = abs(positions_y(j_center_posterior, k_center_posterior)); 
% 
% fprintf('Rest length height summary:\n'); 
% fprintf('Anterior  height = %.3f cm\n', height_anterior)
% fprintf('Posterior height = %.3f cm\n', height_posterior)
% 
% 






