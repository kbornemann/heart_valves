%% BC TUNING

% Parameters
MMHG_TO_CGS = 1333.22368;
cycle_duration = 0.411; 
p_systolic = 60;
p_diastolic = 24;
QPQS = 5; % 4:1

% scale ventricular and atrial data by same factor as aorta
factor_systolic = p_systolic/120;
factor_diastolic = p_diastolic/80;

bump_radius = 0.005; 
n_fourier_coeffs = 600; 
plots = false; 
dt = 2.5e-6;
t_mesh_one_cycle = 0:dt:cycle_duration; 
t = 0:dt:cycle_duration; 

% Rethink
Q_goal_L_per_min = 2 * 0.3504; 
Q_goal_ml_per_s = Q_goal_L_per_min * 1e3 / 60; 
Q_goal_ml_per_cycle = 2 * 2.4; %Q_goal_ml_per_s * cycle_duration; 

ejection_fraction_goal = 0.575;
lv_end_diastolic_volume = 4.2; %Q_goal_ml_per_cycle / ejection_fraction_goal;
lv_volume_initial = lv_end_diastolic_volume - Q_goal_ml_per_cycle/2;
rv_end_diastolic_volume = 4.2;
rv_volume_initial = rv_end_diastolic_volume - Q_goal_ml_per_cycle/2;

% Literature curves: Flowrate mitral
table_q_mitral = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_fig1_q_mitral.csv'); 
times_q_mitral = table_q_mitral.x; 
q_mitral_raw = table_q_mitral.q_mitral; 

% Literature curves: Flowrate tricuspid
table_q_tricuspid = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_fig1_q_mitral.csv'); 
times_q_tricuspid = table_q_tricuspid.x; 
q_tricuspid_raw = table_q_tricuspid.q_mitral; 

% Rescaled time
times_mitral = rescale(times_q_mitral, 0, cycle_duration);
times_tricuspid = rescale(times_q_tricuspid, 0, cycle_duration);

% Flowrate mitral
[a_0_q_mitral, a_n_q_mitral, b_n_q_mitral, Series_q_mitral, ~, ~, Series_q_mitral_derivative] = ...
    series_and_smooth([times_mitral, q_mitral_raw], dt, bump_radius, n_fourier_coeffs, plots); 

vals_series_q_mitral = Series_q_mitral(t);
vals_series_q_mitral_derivative = Series_q_mitral_derivative(t);

% Flowrate tricuspid
[a_0_q_tricuspid, a_n_q_tricuspid, b_n_q_tricuspid, Series_q_tricuspid, ~, ~, Series_q_tricuspid_derivative] = ...
    series_and_smooth([times_tricuspid, q_tricuspid_raw], dt, bump_radius, n_fourier_coeffs, plots); 

vals_series_q_tricuspid = Series_q_tricuspid(t);
vals_series_q_tricuspid_derivative = Series_q_tricuspid_derivative(t);

% Decide whether to process curves from literature (input_literature) or from simulation (input_simulation) 
input_literature = true;
input_simulation = false;

if input_literature

    % Literature curves: pressure 
    table = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_lv_la_ao.csv'); 
    
    times_raw = table.x; 
    pressures_lv_raw = table.pressure_lv; 
    pressures_rv_raw = table.pressure_lv;
    pressures_aorta_raw = table.pressure_aorta;
    pressures_rpa_raw = table.pressure_aorta;
    pressures_lpa_raw = table.pressure_aorta;
    pressures_la_raw = table.pressure_la;
    pressures_ra_raw = table.pressure_la;
    
    % Rescaled time
    times_init = rescale(times_raw, 0, cycle_duration);
    
    % Rescaled LV pressure
    pressures_lv = rescale(pressures_lv_raw, ...
                           factor_diastolic * min(pressures_lv_raw), ...
                           factor_systolic * max(pressures_lv_raw));
    % Rescaled RV pressure
    pressures_rv = rescale(pressures_rv_raw, ...
                           factor_diastolic * min(pressures_rv_raw), ...
                           factor_systolic * max(pressures_rv_raw));

    
    % Rescaled LA pressure
    pressures_la = rescale(pressures_la_raw, ...
                           factor_diastolic * min(pressures_la_raw), ...
                           factor_systolic * max(pressures_la_raw));
    % Rescaled RA pressure
    pressures_ra = rescale(pressures_ra_raw, ...
                           factor_diastolic * min(pressures_ra_raw), ...
                           factor_systolic * max(pressures_ra_raw));
    
    % Rescaled aortic pressure
    pressures_aorta = rescale(pressures_aorta_raw, p_diastolic, p_systolic);

    % Rescaled rpa pressure
    pressures_rpa = rescale(pressures_rpa_raw, p_diastolic, p_systolic);

    % Rescaled lpa pressure
    pressures_lpa = rescale(pressures_lpa_raw, p_diastolic, p_systolic);
    
    rescaled_pressure_plot = true; 
    if rescaled_pressure_plot
        figure; 
        plot(times_init, pressures_lv)
        hold on
        plot(times_init, pressures_rv)
        hold on
        plot(times_init, pressures_la)
        hold on
        plot(times_init, pressures_ra)
        hold on
        plot(times_init, pressures_aorta)
        hold on
        plot(times_init, pressures_rpa)
        hold on
        plot(times_init, pressures_lpa)
        legend('p_lv','p_rv','p_la','p_ra','p_aorta','p_rpa','p_lpa')
    end 
  
    % Define LV pressure and derivative
    [a_0_pressure_lv, a_n_pressure_lv, b_n_pressure_lv, Series_pressure_lv, ~, ~, Series_pressure_lv_derivative] = ... 
        series_and_smooth([times_init, pressures_lv], dt, bump_radius, n_fourier_coeffs, plots); 
    
    vals_series_p_lv = Series_pressure_lv(t); 
    vals_series_p_lv_derivative = Series_pressure_lv_derivative(t); 

    % Define RV pressure and derivative
    [a_0_pressure_rv, a_n_pressure_rv, b_n_pressure_rv, Series_pressure_rv, ~, ~, Series_pressure_rv_derivative] = ... 
        series_and_smooth([times_init, pressures_rv], dt, bump_radius, n_fourier_coeffs, plots); 
    
    vals_series_p_rv = Series_pressure_rv(t); 
    vals_series_p_rv_derivative = Series_pressure_rv_derivative(t); 
    
    % Define aorta pressure and derivative
    [a_0_pressure_aorta, a_n_pressure_aorta, b_n_pressure_aorta, Series_pressure_aorta, ~, ~, Series_pressure_aorta_derivative] = ...
        series_and_smooth([times_init, pressures_aorta], dt, bump_radius, n_fourier_coeffs, plots); 
    
    vals_series_p_aorta = Series_pressure_aorta(t); 
    vals_series_p_aorta_derivative = Series_pressure_aorta_derivative(t); 

    % Define rpa pressure and derivative
    [a_0_pressure_rpa, a_n_pressure_rpa, b_n_pressure_rpa, Series_pressure_rpa, ~, ~, Series_pressure_rpa_derivative] = ...
        series_and_smooth([times_init, pressures_rpa], dt, bump_radius, n_fourier_coeffs, plots); 
    
    vals_series_p_rpa = Series_pressure_rpa(t); 
    vals_series_p_rpa_derivative = Series_pressure_rpa_derivative(t); 

    % Define lpa pressure and derivative
    [a_0_pressure_lpa, a_n_pressure_lpa, b_n_pressure_lpa, Series_pressure_lpa, ~, ~, Series_pressure_lpa_derivative] = ...
        series_and_smooth([times_init, pressures_lpa], dt, bump_radius, n_fourier_coeffs, plots); 
    
    vals_series_p_lpa = Series_pressure_lpa(t); 
    vals_series_p_lpa_derivative = Series_pressure_lpa_derivative(t); 

    % Literature curves: flowrate truncal
    table_q_truncal = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_fig1_q_aorta.csv'); 
    times_q_truncal = table_q_truncal.x; 
    q_truncal_raw = table_q_truncal.q_aorta; 

    % Rescaled time
    times_truncal = rescale(times_q_truncal, 0, cycle_duration);
    
    [a_0_q_truncal, a_n_q_truncal, b_n_q_truncal, Series_q_truncal, ~, ~, Series_q_truncal_derivative] = ...
    series_and_smooth([times_truncal, 2*q_truncal_raw], dt, bump_radius, n_fourier_coeffs, plots); 

    vals_series_q_truncal = Series_q_truncal(t);
    vals_series_q_truncal_derivative = Series_q_truncal_derivative(t);

    % Literature curves: flowrate aorta
    table_q_aorta = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_fig1_q_aorta.csv'); 
    times_q_aorta = table_q_aorta.x; 
    q_aorta_raw = table_q_aorta.q_aorta; 
    q_rpa_raw = table_q_aorta.q_aorta; 
    q_lpa_raw = table_q_aorta.q_aorta; 

    % Rescaled time
    times_aorta = rescale(times_q_aorta, 0, cycle_duration);
    times_rpa = rescale(times_q_aorta, 0, cycle_duration);
    times_lpa = rescale(times_q_aorta, 0, cycle_duration);
    
    [a_0_q_aorta, a_n_q_aorta, b_n_q_aorta, Series_q_aorta, ~, ~, Series_q_aorta_derivative] = ...
    series_and_smooth([times_aorta, q_aorta_raw], dt, bump_radius, n_fourier_coeffs, plots); 

    vals_series_q_aorta = Series_q_aorta(t);
    vals_series_q_aorta_derivative = Series_q_aorta_derivative(t);

    [a_0_q_rpa, a_n_q_rpa, b_n_q_rpa, Series_q_rpa, ~, ~, Series_q_rpa_derivative] = ...
    series_and_smooth([times_rpa, q_rpa_raw], dt, bump_radius, n_fourier_coeffs, plots); 

    vals_series_q_rpa = Series_q_rpa(t);
    vals_series_q_rpa_derivative = Series_q_rpa_derivative(t);

    [a_0_q_lpa, a_n_q_lpa, b_n_q_lpa, Series_q_lpa, ~, ~, Series_q_lpa_derivative] = ...
    series_and_smooth([times_lpa, q_lpa_raw], dt, bump_radius, n_fourier_coeffs, plots); 

    vals_series_q_lpa = Series_q_lpa(t);
    vals_series_q_lpa_derivative = Series_q_lpa_derivative(t);

% elseif input_simulation
% 
%     load('/Users/kmbornemann/KMB_proc/GCF Rev/2-10-25/truncal_preop_opt/BC_tuning/aortic_7217266_preop/bc_data.mat');
% 
%     start = 137000;
% 
%     times = times(start:end);
% 
%     times_sim = rescale(times, 0, cycle_duration);
% 
%     q_lv_sim = q_lvot(start:end);
%     q_rv_sim = q_rvot(start:end);
%     q_truncal_sim = q_lvot(start:end) + q_rvot(start:end);
%     q_aorta_sim = q_aorta(start:end);
%     q_rpa_sim = q_rpa(start:end);
%     q_lpa_sim = q_lpa(start:end);
% 
%     p_lv_sim = p_lvot(start:end);
%     p_rv_sim = p_rvot(start:end);
%     p_aorta_sim = p_aorta(start:end);
%     p_rpa_sim = p_rpa(start:end);
%     p_lpa_sim = p_lpa(start:end);
% 
%     t = times_sim;
% 
%     % LV pressure and derivative
%     [a_0_pressure_lv, a_n_pressure_lv, b_n_pressure_lv, Series_pressure_lv, ~, ~, Series_pressure_lv_derivative] = ... 
%         series_and_smooth([t, p_lv_sim], dt, bump_radius, n_fourier_coeffs, plots); 
% 
%     vals_series_p_lv = Series_pressure_lv(t); 
%     vals_series_p_lv_derivative = Series_pressure_lv_derivative(t); 
% 
% 
%     % RV pressure and derivative
%     [a_0_pressure_rv, a_n_pressure_rv, b_n_pressure_rv, Series_pressure_rv, ~, ~, Series_pressure_rv_derivative] = ... 
%         series_and_smooth([t, p_rv_sim], dt, bump_radius, n_fourier_coeffs, plots); 
% 
%     vals_series_p_rv = Series_pressure_rv(t); 
%     vals_series_p_rv_derivative = Series_pressure_rv_derivative(t); 
% 
%     % Aorta pressure from simulation
%     [a_0_pressure_aorta, a_n_pressure_aorta, b_n_pressure_aorta, Series_pressure_aorta, ~, ~, Series_pressure_aorta_derivative] = ... 
%         series_and_smooth([t, p_aorta_sim], dt, bump_radius, n_fourier_coeffs, plots); 
% 
%     vals_series_p_aorta = Series_pressure_aorta(t); 
%     vals_series_p_aorta_derivative = Series_pressure_aorta_derivative(t); 
% 
%     % RPA pressure from simulation
%     [a_0_pressure_rpa, a_n_pressure_rpa, b_n_pressure_rpa, Series_pressure_rpa, ~, ~, Series_pressure_rpa_derivative] = ... 
%         series_and_smooth([t, p_rpa_sim], dt, bump_radius, n_fourier_coeffs, plots); 
% 
%     vals_series_p_rpa = Series_pressure_rpa(t); 
%     vals_series_p_rpa_derivative = Series_pressure_rpa_derivative(t); 
% 
%     % LPA pressure from simulation
%     [a_0_pressure_lpa, a_n_pressure_lpa, b_n_pressure_lpa, Series_pressure_lpa, ~, ~, Series_pressure_lpa_derivative] = ... 
%         series_and_smooth([t, p_lpa_sim], dt, bump_radius, n_fourier_coeffs, plots); 
% 
%     vals_series_p_lpa = Series_pressure_lpa(t); 
%     vals_series_p_lpa_derivative = Series_pressure_lpa_derivative(t); 
% 
%     % Flowrate truncal from simulation
%     [a_0_q_truncal, a_n_q_truncal, b_n_q_truncal, Series_q_truncal, ~, ~, Series_q_truncal_derivative] = ... 
%         series_and_smooth([t, q_truncal_sim], dt, bump_radius, n_fourier_coeffs, plots); 
% 
%     vals_series_q_truncal = Series_q_truncal(t); 
%     vals_series_q_truncal_derivative = Series_q_truncal_derivative(t); 
% 
%     % Flowrate aorta from simulation
%     [a_0_q_aorta, a_n_q_aorta, b_n_q_aorta, Series_q_aorta, ~, ~, Series_q_aorta_derivative] = ... 
%         series_and_smooth([t, q_aorta_sim], dt, bump_radius, n_fourier_coeffs, plots); 
% 
%     vals_series_q_aorta = Series_q_aorta(t); 
%     vals_series_q_aorta_derivative = Series_q_aorta_derivative(t); 
% 
%     % Flowrate RPA from simulation
%     [a_0_q_rpa, a_n_q_rpa, b_n_q_rpa, Series_q_rpa, ~, ~, Series_q_rpa_derivative] = ... 
%         series_and_smooth([t, q_rpa_sim], dt, bump_radius, n_fourier_coeffs, plots); 
% 
%     vals_series_q_rpa = Series_q_rpa(t); 
%     vals_series_q_rpa_derivative = Series_q_rpa_derivative(t); 
% 
%     % Flowrate LPA from simulation
%     [a_0_q_lpa, a_n_q_lpa, b_n_q_lpa, Series_q_lpa, ~, ~, Series_q_lpa_derivative] = ... 
%         series_and_smooth([t, q_lpa_sim], dt, bump_radius, n_fourier_coeffs, plots); 
% 
%     vals_series_q_lpa = Series_q_lpa(t); 
%     vals_series_q_lpa_derivative = Series_q_lpa_derivative(t); 
% 
%     rescaled_pressure_plot = true; 
%     if rescaled_pressure_plot
%         figure; 
%         plot(t', vals_series_p_lv)
%         hold on
%         plot(t', vals_series_p_rv)
%         hold on
%         plot(t', vals_series_p_aorta)
%         hold on
%         plot(t', vals_series_p_rpa)
%         hold on
%         plot(t', vals_series_p_lpa)
%         legend('p_lv','p_rv','p_aorta','p_rpa','p_rpa')
%     end 
% 

end

% Normalize flow rate integrals to be equal and to desired target

q_truncal_cumulative = dt * trapz(vals_series_q_truncal)
q_mitral_cumulative = dt * trapz(vals_series_q_mitral) 
q_tricuspid_cumulative = dt * trapz(vals_series_q_tricuspid)

q_aorta_cumulative = dt * trapz(vals_series_q_aorta)
q_rpa_cumulative = dt * trapz(vals_series_q_rpa)
q_lpa_cumulative = dt * trapz(vals_series_q_lpa)

scaling_q_mitral = (Q_goal_ml_per_cycle/2) / q_mitral_cumulative;
scaling_q_tricuspid = (Q_goal_ml_per_cycle/2) / q_tricuspid_cumulative;
scaling_q_truncal = Q_goal_ml_per_cycle / q_truncal_cumulative;
scaling_q_aorta = ((1/QPQS)*Q_goal_ml_per_cycle) / q_aorta_cumulative;
scaling_q_rpa = ((2/QPQS)*Q_goal_ml_per_cycle) / q_rpa_cumulative;
scaling_q_lpa = ((2/QPQS)*Q_goal_ml_per_cycle) / q_lpa_cumulative;

Series_q_mitral_scaled = @(t) scaling_q_mitral * Series_q_mitral(t);
vals_series_q_mitral_scaled = Series_q_mitral_scaled(t);

Series_q_tricuspid_scaled = @(t) scaling_q_tricuspid * Series_q_tricuspid(t);
vals_series_q_tricuspid_scaled = Series_q_tricuspid_scaled(t);

Series_q_truncal_scaled = @(t) scaling_q_truncal * Series_q_truncal(t);
vals_series_q_truncal_scaled = Series_q_truncal_scaled(t);

Series_q_aorta_scaled = @(t) scaling_q_aorta * Series_q_aorta(t);
vals_series_q_aorta_scaled = Series_q_aorta_scaled(t);

Series_q_rpa_scaled = @(t) scaling_q_rpa * Series_q_rpa(t);
vals_series_q_rpa_scaled = Series_q_rpa_scaled(t);

Series_q_lpa_scaled = @(t) scaling_q_lpa * Series_q_lpa(t);
vals_series_q_lpa_scaled = Series_q_lpa_scaled(t);

q_mitral_cumulative_scaled = dt * trapz(vals_series_q_mitral_scaled)
q_tricuspid_cumulative_scaled = dt * trapz(vals_series_q_tricuspid_scaled)
q_truncal_cumulative_scaled = dt * trapz(vals_series_q_truncal_scaled)
q_aorta_cumulative_scaled = dt * trapz(vals_series_q_aorta_scaled)
q_rpa_cumulative_scaled = dt * trapz(vals_series_q_rpa_scaled)
q_lpa_cumulative_scaled = dt * trapz(vals_series_q_lpa_scaled)

% OUTPUT FILE: FOURIER_COEFFS_Q_IN - fourier_coeffs_Q_mi.txt
base_name = 'fourier_coeffs';
suffix = '_Q_mi';
file_name = strcat(base_name, suffix, '.txt'); 
output_series_coeffs_to_txt(scaling_q_mitral * a_0_q_mitral, scaling_q_mitral * a_n_q_mitral, scaling_q_mitral * b_n_q_mitral, n_fourier_coeffs, cycle_duration, file_name); 

% OUTPUT FILE: FOURIER_COEFFS_Q_IN - fourier_coeffs_Q_tri.txt
base_name = 'fourier_coeffs';
suffix = '_Q_tri';
file_name = strcat(base_name, suffix, '.txt'); 
output_series_coeffs_to_txt(scaling_q_tricuspid * a_0_q_tricuspid, scaling_q_tricuspid * a_n_q_tricuspid, scaling_q_tricuspid * b_n_q_tricuspid, n_fourier_coeffs, cycle_duration, file_name); 

% Convert pressure units to CGS
Series_pressure_lv_cgs = @(t) MMHG_TO_CGS * Series_pressure_lv(t);
Series_pressure_lv_derivative_cgs = @(t) MMHG_TO_CGS * Series_pressure_lv_derivative(t);
vals_series_pressure_lv_cgs = Series_pressure_lv_cgs(t);
vals_series_pressure_lv_derivative_cgs = Series_pressure_lv_derivative_cgs(t);

Series_pressure_rv_cgs = @(t) MMHG_TO_CGS * Series_pressure_rv(t);
Series_pressure_rv_derivative_cgs = @(t) MMHG_TO_CGS * Series_pressure_rv_derivative(t);
vals_series_pressure_rv_cgs = Series_pressure_rv_cgs(t);
vals_series_pressure_rv_derivative_cgs = Series_pressure_rv_derivative_cgs(t);

Series_pressure_aorta_cgs = @(t) MMHG_TO_CGS * Series_pressure_aorta(t);
Series_pressure_aorta_derivative_cgs = @(t) MMHG_TO_CGS * Series_pressure_aorta_derivative(t);
vals_series_pressure_aorta_cgs = Series_pressure_aorta_cgs(t);
vals_series_pressure_aorta_derivative_cgs = Series_pressure_aorta_derivative_cgs(t);

Series_pressure_rpa_cgs = @(t) MMHG_TO_CGS * Series_pressure_rpa(t);
Series_pressure_rpa_derivative_cgs = @(t) MMHG_TO_CGS * Series_pressure_rpa_derivative(t);
vals_series_pressure_rpa_cgs = Series_pressure_rpa_cgs(t);
vals_series_pressure_rpa_derivative_cgs = Series_pressure_rpa_derivative_cgs(t);

Series_pressure_lpa_cgs = @(t) MMHG_TO_CGS * Series_pressure_lpa(t);
Series_pressure_lpa_derivative_cgs = @(t) MMHG_TO_CGS * Series_pressure_lpa_derivative(t);
vals_series_pressure_lpa_cgs = Series_pressure_lpa_cgs(t);
vals_series_pressure_lpa_derivative_cgs = Series_pressure_lpa_derivative_cgs(t);

% LV volume
vals_lv_volume = zeros(size(vals_series_q_truncal)); 
vals_lv_volume_deriv = zeros(size(vals_series_q_truncal)); 
vals_lv_volume(1) = lv_volume_initial;

% RV volume
vals_rv_volume = zeros(size(vals_series_q_truncal)); 
vals_rv_volume_deriv = zeros(size(vals_series_q_truncal)); 
vals_rv_volume(1) = rv_volume_initial;

% LV
for j=2:length(vals_lv_volume)
    vals_lv_volume(j) = vals_lv_volume(j-1) + dt * (vals_series_q_mitral_scaled(j) - vals_series_q_truncal_scaled(j)); 
end 

for j=1:length(vals_rv_volume)
    vals_lv_volume_deriv(j) = vals_series_q_mitral_scaled(j) - vals_series_q_truncal_scaled(j); 
end 

% RV
for j=2:length(vals_rv_volume)
    vals_rv_volume(j) = vals_rv_volume(j-1) + dt * (vals_series_q_tricuspid_scaled(j) - vals_series_q_truncal_scaled(j)); 
end 

for j=1:length(vals_rv_volume)
    vals_rv_volume_deriv(j) = vals_series_q_tricuspid_scaled(j) - vals_series_q_truncal_scaled(j); 
end 

two_hill = true; 
if two_hill

    % Insert parameters from svzerod_lv_tune/run_lv_tune.py
    % Brown AMBE 2023 

    t_shift =  0.050690005643485476
    tau_1 =  0.18939919642752506
    tau_2 =  0.3542802482730453
    m1 =  18.25804089495962
    m2 =  31.316768151285913

    g1 = @(t) (t./tau_1).^m1; 
    g2 = @(t) (t./tau_2).^m2; 
    
    r1 = @(t) g1(t) ./ (1 + g1(t)); 
    r2 = @(t) 1 ./ (1 + g2(t)); 
    two_hill_product = @(t) (g1(t) ./ (1 + g1(t))) .* (1 ./ (1 + g2(t))); 

    times_hill = 0:dt:cycle_duration; 
    two_hill_product_maximum = max(two_hill_product(times_hill)); 
    
    k_coeff_two_hill = 1 / two_hill_product_maximum; 

    t_shift_hill = t_shift;

    two_hill_function = @(t) k_coeff_two_hill * two_hill_product(mod(t - t_shift_hill,cycle_duration)); 

    figure; 
    hold on    
    plot(times_hill, two_hill_function(times_hill)); 
    plot(times_hill, r1(times_hill))
    plot(times_hill, r2(times_hill))    
    legend('hill function', 'r1', 'r2')
    title('hill function')
    
end 

% Two-hill function
[a_0_activation_two_hill, a_n_activation_two_hill, b_n_activation_two_hill, Series_activation_two_hill] = ... 
    series_and_smooth([t', two_hill_function(t)'], dt, bump_radius, n_fourier_coeffs, plots); 

vals_series_activation_two_hill = Series_activation_two_hill(t);

% OUTPUT FILE: FOURIER_COEFFS_ACT - fourier_coeffs_lv_activation_two_hill.txt
base_name = 'fourier_coeffs';
suffix = '_lv_activation_two_hill'
file_name = strcat(base_name, suffix, '.txt'); 
output_series_coeffs_to_txt(a_0_activation_two_hill, a_n_activation_two_hill, b_n_activation_two_hill, n_fourier_coeffs, cycle_duration, file_name); 

% % OUTPUT FILE: FOURIER_COEFFS_ACT - fourier_coeffs_rv_activation_two_hill.txt
% base_name = 'fourier_coeffs';
% suffix = '_rv_activation_two_hill'
% file_name = strcat(base_name, suffix, '.txt'); 
% output_series_coeffs_to_txt(a_0_activation_two_hill, a_n_activation_two_hill, b_n_activation_two_hill, n_fourier_coeffs, cycle_duration, file_name); 

% OUTPUT FILE: NOT USED - fourier_coeffs_lv_activation.txt
p_lv_activation_threshold = 10; 
activation_data_unscaled = (vals_series_p_lv' > p_lv_activation_threshold) .* vals_series_p_lv';
activation_data = activation_data_unscaled / max(activation_data_unscaled);

[a_0_activation, a_n_activation, b_n_activation, Series_activation] = ... 
    series_and_smooth([t', activation_data], dt, bump_radius, n_fourier_coeffs, plots); 

vals_series_activation = Series_activation(t);

base_name = 'fourier_coeffs';
suffix = '_lv_activation'
file_name = strcat(base_name, suffix, '.txt'); 
output_series_coeffs_to_txt(a_0_activation, a_n_activation, b_n_activation, n_fourier_coeffs, cycle_duration, file_name); 

% base_name = 'fourier_coeffs';
% suffix = '_rv_activation'
% file_name = strcat(base_name, suffix, '.txt'); 
% output_series_coeffs_to_txt(a_0_activation, a_n_activation, b_n_activation, n_fourier_coeffs, cycle_duration, file_name); 


series_plots = true; 
if series_plots
    figure; 
    plot(t, vals_series_p_lv)
    hold on
    plot(t, vals_series_p_rv)
    hold on
    plot(t, vals_series_p_aorta)
    hold on
    plot(t, vals_series_p_rpa)
    hold on
    plot(t, vals_series_p_lpa)
    legend('p_lv', 'p_rv', 'p_aorta', 'p_rpa', 'p_lpa')
    title('Lv pressure, RV pressure, aorta pressure, rpa pressure, lpa pressure')

    figure;
    plot(t, vals_series_q_mitral_scaled)
    hold on
    plot(t, vals_series_q_tricuspid_scaled)
    hold on
    plot(t, vals_series_q_aorta_scaled)
    hold on
    plot(t, vals_series_q_rpa_scaled)
    hold on
    plot(t, vals_series_q_lpa_scaled)
    legend('q_mitral', 'q_tricuspid', 'q_aorta', 'q_rpa', 'q_lpa')
    title('q_mitral, q_tricuspid, q_aorta, q_rpa, q_lpa')

    figure;
    plot(t, vals_lv_volume)
    hold on
    plot(t, vals_rv_volume)
    title('lv volume, rv volume')

    %figure;
    %plot(t, vals_series_activation)
    %title('vals_series_activation')
end 



output_to_sv0d = true; 
if output_to_sv0d
    format long 
    f = fopen("array_values.txt", "w");
    fprintf(f, '        "y": {\n');
    print_var_string(f,t,'flow:INFLOW_LV:vessel_LV', 0.5.*vals_series_q_truncal_scaled)
    print_var_string(f,t,'pressure:INFLOW_LV:vessel_LV', vals_series_pressure_lv_cgs)
    print_var_string(f,t,'flow:INFLOW_RV:vessel_RV', 0.5.*vals_series_q_tricuspid_scaled)
    print_var_string(f,t,'pressure:INFLOW_RV:vessel_RV', vals_series_pressure_rv_cgs)
    print_var_string(f,t,'flow:vessel_aorta:OUTLET_AORTA', vals_series_q_aorta_scaled)
    print_var_string(f,t,'pressure:vessel_aorta:OUTLET_AORTA', vals_series_pressure_aorta_cgs)
    print_var_string(f,t,'flow:vessel_rpa:OUTLET_RPA', vals_series_q_rpa_scaled)
    print_var_string(f,t,'pressure:vessel_rpa:OUTLET_RPA', (1/2.25).*vals_series_pressure_aorta_cgs)
    print_var_string(f,t,'flow:vessel_lpa:OUTLET_LPA', vals_series_q_lpa_scaled)
    print_var_string(f,t,'pressure:vessel_lpa:OUTLET_LPA', (1/2.25).*vals_series_pressure_aorta_cgs)
    print_var_string(f,t,'t', t');
end 

function print_var_string(f,t,name,vals)

    fprintf(f, '        "%s": [\n', name);
    fprintf(f, '            ');
    for j = 1:length(t)
        fprintf(f, '%.14f', vals(j));
        if j < length(t)
            fprintf(f, ', ');
        end 
    end 
    fprintf(f, '\n');
    fprintf(f, '        ],\n');

end 

% % run the lpn 
% dt_lpn = 1e-4; 
% n_cycles = 1; 
% t_final = cycle_duration * n_cycles;
% P_ao_initial = 94*MMHG_TO_CGS;
% R_proximal = 83.6698220729; 
% C =  0.00167055364456;
% R_distal = 1287.64596307;
% 
% v_initial = ventricular_volume_initial;
% 
% KERCKHOFFS = true; 
% 
% if KERCKHOFFS 
%     Vrd = 26.1; % ml 
%     Vrs = 18; % ml
% 
%     % parameters from KERCKHOFFS AMBE 2006 
%     % lv 
%     % note that min compliance is actually the larger value under confusing naming convention 
%     C_min_ml_over_kPa = 11.0; 
%     C_max_ml_over_kPa = 0.946; 
% 
%     % E_min_kPa = 1/C_min_ml_over_kPa; 
%     % E_max_kPa = 1/C_max_ml_over_kPa; 
% 
%     ML_OVER_KPA_TO_ML_OVER_DYNEPERCM2 = 1e-4; 
% 
%     C_min_scaling = 5; 
%     C_max_scaling = 5; 
% 
%     C_min_ml_over_dynespercm2 = C_min_scaling * C_min_ml_over_kPa * ML_OVER_KPA_TO_ML_OVER_DYNEPERCM2; 
%     C_max_ml_over_dynespercm2 = C_max_scaling * C_max_ml_over_kPa * ML_OVER_KPA_TO_ML_OVER_DYNEPERCM2; 
% 
%     Emax = 1/C_max_ml_over_dynespercm2
%     Emin = 1/C_min_ml_over_dynespercm2
% 
%     Emin_mmHg_over_ml = Emin / MMHG_TO_CGS 
%     Emax_mmHg_over_ml = Emax / MMHG_TO_CGS
% 
% else 
%     error('not implemented');
% end 
% 
% R_av_closed = 100000;
% steepness_av = 0.00001; 
% 
% [times_lpn, P_lv, Q_ao, P_ao, V_lv, R_tanh] = solve_lv_ao_lpn(dt_lpn, t_final, v_initial, Vrd, Vrs, Emax, Emin, ...
%                                      Series_q_mitral_scaled, Series_activation, ...
%                                      P_ao_initial, R_proximal, C, R_distal, r_av, R_av_closed, steepness_av);
% 
% figure;
% subplot(6,1,1)
% plot(t, vals_series_pressure_lv)
% hold on 
% plot(t, vals_series_pressure_aorta)
% plot(times_lpn, P_lv/MMHG_TO_CGS); 
% plot(times_lpn, P_ao/MMHG_TO_CGS);
% xlim([0 max(times_lpn)])
% 
% legend('lv exp', 'ao exp', 'lv lpn', 'ao lpn');
% 
% subplot(6,1,2)
% hold on 
% plot(t, vals_series_q_mitral_scaled)
% plot(t, vals_series_q_aorta_scaled)
% plot(times_lpn, Q_ao)
% xlim([0 max(times_lpn)])
% legend('q mi exp', 'q ao exp', 'q ao lpn')
% 
% subplot(6,1,3)
% plot(times_hill, two_hill_function(times_hill)); 
% 
% subplot(6,1,4)
% plot(t, vals_ventricular_volume)
% hold on 
% plot(times_lpn, V_lv)
% xlim([0 max(times_lpn)])
% legend('V integrated', 'V lpn')
% 
% subplot(6,1,5);
% plot(t, vals_series_activation_two_hill);
% hold on 
% plot(t, vals_series_activation);
% xlim([0 max(times_lpn)])
% legend('two hill', 'act from p lv')
% title('act')
% 
% subplot(6,1,6)
% plot(times_lpn, R_tanh)
% xlim([0 max(times_lpn)])
% title('r valve')
% figure; 
% hold on 
% plot(V_lv, P_lv/MMHG_TO_CGS)
% xlabel('V ml')
% ylabel('P LV mmHg')









