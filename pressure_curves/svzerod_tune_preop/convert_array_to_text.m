
format long 
f = fopen("array_values.txt", "w");
t = times_lim;
fprintf(f, '        "y": {\n');
print_var_string(f,t,'p_lvot', p_lvot_lim)
print_var_string(f,t,'p_rvot', p_rvot_lim)
print_var_string(f,t,'times', times_lim)


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