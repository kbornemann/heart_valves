function F = difference_equations_aortic(leaflet)
    % 
    % Evaluation of the global difference equations at j,k
    % Requires reference configuration R 
    % 
    % Input
    %     leaflet    Current parameters 
    %
    % Output
    %     F        Values of all difference equation, 3 by triangular array 
    % 

    % Copyright (c) 2019, Alexander D. Kaiser
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are met:
    % 
    % 1. Redistributions of source code must retain the above copyright notice, this
    %    list of conditions and the following disclaimer.
    % 
    % 2. Redistributions in binary form must reproduce the above copyright notice,
    %    this list of conditions and the following disclaimer in the documentation
    %    and/or other materials provided with the distribution.
    % 
    % 3. Neither the name of the copyright holder nor the names of its
    %    contributors may be used to endorse or promote products derived from
    %    this software without specific prior written permission.
    % 
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    % DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
    % FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    % DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    % SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    % CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    % OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    % OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    X_current              = leaflet.X; 
    p_0                    = leaflet.p_0; 
    alpha                  = leaflet.alpha; 
    beta                   = leaflet.beta; 
    c_dec_radial           = leaflet.c_dec_radial; 
    c_dec_circumferential  = leaflet.c_dec_circumferential; 
    j_max                  = leaflet.j_max; 
    k_max                  = leaflet.k_max; 
    du                     = leaflet.du; 
    is_internal            = leaflet.is_internal; 
    is_bc                  = leaflet.is_bc; 
    
    if isfield(leaflet, 'periodic_j')
        periodic_j = leaflet.periodic_j; 
    else
        periodic_j = zeros(k_max,1); 
    end 
    
    if isfield(leaflet, 'decreasing_tension') && leaflet.decreasing_tension
        decreasing_tension    = true;  
    else 
        decreasing_tension    = false;  
    end 
    

    if isfield(leaflet, 'tension_debug') && leaflet.tension_debug
        tension_debug = true; 
    else 
        tension_debug = false; 
    end 

    
    F_leaflet = zeros(size(X_current)); 
    
    % Internal leaflet part 
    for j=1:j_max
        for k=1:k_max
            if is_internal(j,k) 

                X = X_current(:,j,k); 

                F_tmp = zeros(3,1);

                % pressure term first  
                if p_0 ~= 0
                    
                    [j_plus__1 j_minus_1 k_plus__1 k_minus_1 m] = get_pressure_nbrs(leaflet,j,k); 
                    
                    F_tmp = F_tmp + p_0 * m * cross(X_current(:,j_plus__1,k) - X_current(:,j_minus_1,k), X_current(:,j,k_plus__1) - X_current(:,j,k_minus_1));                     
                end 

                
                % u type fibers 
                for j_nbr_tmp = [j-1,j+1]
                    
                    k_nbr_tmp = k; 
                    
                    [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                    
                    if valid

                        X_nbr = X_current(:,j_nbr,k_nbr); 

                        alpha_tmp     = alpha(j_spr,k_spr); 
                        c_dec_tension = c_dec_circumferential(j_spr,k_spr); 
                        
                        tension = alpha_tmp;  

                        if decreasing_tension && (alpha_tmp ~= 0)
                            tension = tension + alpha_tmp * tension_decreasing(X, X_nbr, du, c_dec_tension) ; 
                        end 

                        if tension_debug
                            dec = tension_decreasing(X, X_nbr, du, c_dec_tension) ; 
                            fprintf('tension = %e, dec_tension = %f, (j,k) = (%d, %d) circ\n', tension, dec, j, k); 
                        end 

                        F_tmp = F_tmp + du * tension * (X_nbr-X)/norm(X_nbr-X); 
                    
                    end 
                    
                    
                    
                end 

                % v type fibers 
                for k_nbr_tmp = [k-1,k+1]

                    j_nbr_tmp = j; 
                    
                    [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                    
                    if valid
                    
                        X_nbr = X_current(:,j_nbr,k_nbr); 

                        beta_tmp      = beta(j_spr,k_spr); 
                        c_dec_tension = c_dec_radial(j_spr,k_spr); 
                        
                        tension = beta_tmp; 

                        if decreasing_tension && (beta_tmp ~= 0)
                            tension = tension + beta_tmp * tension_decreasing(X, X_nbr, du, c_dec_tension) ; 
                        end

                        if tension_debug
                            dec = tension_decreasing(X, X_nbr, du, c_dec_tension) ; 
                            fprintf('tension = %e, dec_tension = %f, (j,k) = (%d, %d) radial\n', tension, dec, j, k); 
                        end 
                        
                        F_tmp = F_tmp + du * tension * (X_nbr-X)/norm(X_nbr-X); 
                                               
                    end 

                end 
                
                F_leaflet(:,j,k) = F_tmp;

            end
        end 
    end
    
    F = linearize_internal_points(leaflet, F_leaflet); 
    
end 



