run("bc_data.m")

MMHG_TO_CGS = 1333.22368;

% table_q_mitral = readtable('pressure_curve_data/Physiologic_Mechanisms_Aortic_Insufficiency_Yellin/PhysiologicMechanismsinAorticInsufficiency_fig1_q_mitral.csv'); 
% times_q_mitral = table_q_mitral.x; 
% q_mitral_raw = table_q_mitral.q_mitral; 

dt = times(2)-times(1);
cycle_duration = 0.382;
bump_radius = 0.01;
n_fourier_coeffs = 600;
plots=false;
t = 0:dt:cycle_duration;

times = times(1:152801);
q_aorta = q_aorta(1:152801);
q_lv = q_ventricle(1:152801);
p_aorta = p_aorta(1:152801);
p_lv = p_lv(1:152801);

% % Flow aorta
% [a_0_q_aorta, a_n_q_aorta, b_n_q_aorta, Series_q_aorta, ~, ~, Series_q_aorta_derivative] = ...
%     series_and_smooth([times, q_aorta], dt, bump_radius, n_fourier_coeffs, plots); 
% 
% 
% vals_series_q_aorta = Series_q_aorta(t);
% vals_series_q_aorta_derivative = Series_q_aorta_derivative(t);
% 
% % Flow LV
% [a_0_q_lv, a_n_q_lv, b_n_q_lv, Series_q_lv, ~, ~, Series_q_lv_derivative] = ...
%     series_and_smooth([times, q_lv], dt, bump_radius, n_fourier_coeffs, plots); 
% 
% 
% vals_series_q_lv = Series_q_lv(t);
% vals_series_q_lv_derivative = Series_q_lv_derivative(t);
% 
% 
% % Pressure LV
% [a_0_p_lv, a_n_p_lv, b_n_p_lv, Series_p_lv, ~, ~, Series_p_lv_derivative] = ...
%     series_and_smooth([times, p_lv], dt, bump_radius, n_fourier_coeffs, plots); 
% 
% Series_p_lv_cgs = @(t) MMHG_TO_CGS * Series_p_lv(t);
% Series_p_lv_derivative_cgs = @(t) MMHG_TO_CGS * Series_p_lv_derivative(t);
% vals_series_p_lv_cgs = Series_p_lv_cgs(t); 
% vals_series_p_lv_derivative_cgs = Series_p_lv_derivative_cgs(t); 
% 
% % Pressure aorta
% [a_0_p_aorta, a_n_p_aorta, b_n_p_aorta, Series_p_aorta, ~, ~, Series_p_aorta_derivative] = ...
%     series_and_smooth([times, p_aorta], dt, bump_radius, n_fourier_coeffs, plots); 
% 
% Series_p_aorta_cgs = @(t) MMHG_TO_CGS * Series_p_aorta(t);
% Series_p_aorta_derivative_cgs = @(t) MMHG_TO_CGS * Series_p_aorta_derivative(t);
% vals_series_p_aorta_cgs = Series_p_aorta_cgs(t); 
% vals_series_p_aorta_derivative_cgs = Series_p_aorta_derivative_cgs(t); 


output_to_sv0d = true; 
if output_to_sv0d
    % output 
    format long 

    f = fopen("array_values.txt", "w");

    print_var_string(f,t,'flow:ventricle:valve1', q_aorta); %vals_series_q_aorta)
    print_var_string(f,t,'pressure:ventricle:valve1', MMHG_TO_CGS * p_lv); %vals_series_p_lv_cgs)
    print_var_string(f,t,'flow:vessel:OUTLET', q_aorta); %vals_series_q_aorta)
    print_var_string(f,t,'pressure:vessel:OUTLET', MMHG_TO_CGS * p_aorta); %vals_series_p_aorta_cgs)
    print_var_string(f,t,'Q', q_lv);
    print_var_string(f,t,'t', times);
    %print_var_string(f,t,'Vc:ventricle', vals_ventricular_volume);

    fprintf(f, '    },\n');

    % fprintf(f, '    "dy": {\n');
    % print_var_string(f,t,'flow:ventricle:valve1', vals_series_q_aorta_derivative)
    % print_var_string(f,t,'pressure:ventricle:valve1', vals_series_p_lv_derivative_cgs)
    % print_var_string(f,t,'flow:vessel:OUTLET', vals_series_q_aorta_derivative)
    % print_var_string(f,t,'pressure:vessel:OUTLET', vals_series_p_aorta_derivative_cgs)
    % %print_var_string(f,t,'Vc:ventricle', vals_ventricular_volume_deriv);
    % fprintf(f, '    },\n');

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