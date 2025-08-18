% Taken from beat 1, p. 227, fig 2 
% 'dynamics of left ventricular filling' Edward Yellin 
% In Cardiac Mechanics and Function in the Normal and Diseased Heart 

cycle_length = 0.411;
p_systolic = 60;
p_diastolic = 24;

% scale ventricular and atrial data by same factor as aorta
factor_systolic = p_systolic/120;
factor_diastolic = p_diastolic/80;


% quadrature spacing 
debug = false; 
if debug 
    dt = 5e-5; 
else 
    dt = 5e-6; 
end 

plots = true; 

base_name = 'fourier_coeffs';

points_one_cycle_lvot = [0.0,   0;
0.02, 0;
0.06, 2; 
0.40, 6; 
0.53, 14; 
0.58, 120; 
0.75, 130; 
0.8, 8]; 

points_one_cycle_rvot = [0.0,   0;
0.02, 0;
0.06, 2; 
0.40, 6; 
0.53, 14; 
0.58, 120; 
0.75, 130; 
0.8, 8]; 

points_one_cycle_atrium = [0.0, 24.555; 
0.06, 4; 
0.40, 7; 
0.47, 20; 
0.53, 5; 
0.58, 7; 
0.7,  10; 
0.8, 24.555]; 

points_one_cycle_aorta = [0.0,  120;  
0.569, 80; 
0.58, 80; 
0.59, 98; 
0.60, 106;
0.61, 115;
0.72, 124;
0.73, 122;
0.745, 119; 
0.755, 100; 
0.76, 122;
0.8, 120]; 

points_one_cycle_lvot(:,1) = rescale(points_one_cycle_lvot(:,1), 0, cycle_length);
points_one_cycle_lvot(:,2) = rescale(points_one_cycle_lvot(:,2), ...
                                          factor_diastolic * min(points_one_cycle_lvot(:,2)), ...
                                          factor_systolic * max(points_one_cycle_lvot(:,2)));

points_one_cycle_rvot(:,1) = rescale(points_one_cycle_rvot(:,1), 0, cycle_length);
points_one_cycle_rvot(:,2) = rescale(points_one_cycle_rvot(:,2), ...
                                          factor_diastolic * min(points_one_cycle_rvot(:,2)), ...
                                          factor_systolic * max(points_one_cycle_rvot(:,2)));

points_one_cycle_atrium(:,1) = rescale(points_one_cycle_atrium(:,1), 0, cycle_length);
points_one_cycle_atrium(:,2) = rescale(points_one_cycle_atrium(:,2), ...
                                          factor_diastolic * min(points_one_cycle_atrium(:,2)), ...
                                          factor_systolic * max(points_one_cycle_atrium(:,2)));

points_one_cycle_aorta(:,1) = rescale(points_one_cycle_aorta(:,1), 0, cycle_length);
points_one_cycle_aorta(:,2) = rescale(points_one_cycle_aorta(:,2), p_diastolic, p_systolic);

suffix = ''; 

file_name = strcat(base_name, suffix, '.txt'); 

bump_radius = .025; 
n_fourier_coeffs = 600; 
% plots = false; 

[a_0_lvot a_n_lvot b_n_lvot Series_lvot] = series_and_smooth(points_one_cycle_lvot, dt, bump_radius, n_fourier_coeffs, plots); 
[a_0_rvot a_n_rvot b_n_rvot Series_rvot] = series_and_smooth(points_one_cycle_rvot, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_lvot_series = Series_lvot(t); 
vals_rvot_series = Series_rvot(t); 

[a_0_atrium a_n_atrium b_n_atrium Series_atrium] = series_and_smooth(points_one_cycle_atrium, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_atrium_series = Series_atrium(t); 

bump_radius = .01; 

[a_0_aorta a_n_aorta b_n_aorta Series_aorta] = series_and_smooth(points_one_cycle_aorta, dt, bump_radius, n_fourier_coeffs, plots); 

t = 0:dt:cycle_length; 
vals_aorta_series = Series_aorta(t); 

min_pressure_aorta = min(vals_aorta_series)
max_pressure_aorta = max(vals_aorta_series)
mean_pressure_aorta = mean(vals_aorta_series)

time_zero_pressure_aorta = vals_aorta_series(1)

vals_plus_one  = [vals_aorta_series(2:end), vals_aorta_series(1)]; 
vals_minus_one = [vals_aorta_series(end), vals_aorta_series(1:(end-1))]; 
dp_dt = (vals_plus_one - vals_minus_one)/(2*dt); 

min_dp_dt_aorta = min(dp_dt)
max_dp_dt_aorta = max(dp_dt)

% fig = figure; 
% plot(t, dp_dt, 'k'); 
% hold on
% title('dP/dt')
% xlabel('t')
% ylabel('dp/dt (mmHg/s)')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% saveas(fig, strcat('dp_dt', suffix))
% 
% fig = figure; 
% plot(t, vals_ventricle_series, 'k'); 
% hold on
% plot(t, vals_atrium_series, '--k'); 
% plot(t, vals_aorta_series, ':k'); 
% title('Three pressures')
% xlabel('t')
% ylabel('p (mmHg)')
% legend('LV', 'LA', 'Aorta','Location', 'SouthEast')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% saveas(fig, strcat('both_pressure_yellin', suffix))
% 
% 
% fig = figure; 
% plot(t, vals_ventricle_series - vals_aorta_series, 'k'); 
% hold on
% plot(t, 0*vals_ventricle_series, 'k'); 
% xlabel('t')
% ylabel('p (mmHg)')
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% title('LV minus aorta')
% saveas(fig, 'pressure_difference_lv_aorta')
% 
% 
% t = 0:dt:(3*cycle_length); 
% vals_ventricle_series = Series_ventricle(t); 
% fig = figure; 
% plot(t, vals_ventricle_series, 'k'); 
% title('Pressures')
% xlabel('t')
% ylabel('p (mmHg)')
% hold on 
% vals_atrium_series = Series_atrium(t); 
% plot(t, vals_atrium_series, '--k');
% vals_aorta_series = Series_aorta(t); 
% plot(t, vals_aorta_series, ':k'); 
% %axis([0 2.4 -10 140])
% set(fig, 'Position', [100, 100, 1000, 500])
% set(fig,'PaperPositionMode','auto')
% legend('Left Ventricle', 'Left Atrium', 'Location', 'NorthWest')
% saveas(fig, 'both_pressure_yellin_three_cycles') 


n_coeffs_to_output = n_fourier_coeffs; 

output_series_coeffs_to_txt(a_0_lvot, a_n_lvot, b_n_lvot, n_coeffs_to_output, cycle_length, 'fourier_coeffs_lvot.txt'); 
output_series_coeffs_to_txt(a_0_rvot, a_n_rvot, b_n_rvot, n_coeffs_to_output, cycle_length, 'fourier_coeffs_rvot.txt'); 
output_series_coeffs_to_txt(a_0_atrium,    a_n_atrium,    b_n_atrium,    n_coeffs_to_output, cycle_length, 'fourier_coeffs_atrium.txt'); 
output_series_coeffs_to_txt(a_0_aorta,     a_n_aorta,     b_n_aorta,     n_coeffs_to_output, cycle_length, 'fourier_coeffs_aorta.txt'); 


% a_0 = a_0_aorta - a_0_ventricle; 
% a_n = a_n_aorta - a_n_ventricle; 
% b_n = b_n_aorta - b_n_ventricle; 
% 
% output_series_coeffs_to_txt(a_0, a_n, b_n, n_coeffs_to_output, cycle_length, 'fourier_coeffs_aorta_minus_ventricle.txt'); 
% 


% close all 
% save(strcat('series_data', suffix)); 






