function F = difference_equations_with_reference(leaflet)
    % 
    % Evaluation of the global difference equations at j,k
    % Uses linear constitutive laws 
    % Requires rest lengths and spring constants to be previously
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

    X_current          = leaflet.X; 
    p_0                = leaflet.p_0; 
    chordae            = leaflet.chordae;
    chordae_idx        = leaflet.chordae_idx; 
    j_max              = leaflet.j_max; 
    k_max              = leaflet.k_max; 
    is_internal        = leaflet.is_internal; 
    is_bc              = leaflet.is_bc; 
    num_trees          = leaflet.num_trees; 
    R_u                = leaflet.R_u;
    k_u                = leaflet.k_u;
    R_v                = leaflet.R_v;
    k_v                = leaflet.k_v;
    
    
    if isfield(leaflet, 'periodic_j')
        periodic_j = leaflet.periodic_j; 
    else
        periodic_j = zeros(k_max,1); 
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

                        tension = tension_with_reference(X, X_nbr, R_u(j_spr,k_spr), k_u(j_spr,k_spr), leaflet); 
                        F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X);                     
                    end 
                    
                end 

                % v type fibers 
                for k_nbr_tmp = [k-1,k+1]

                    j_nbr_tmp = j; 
                    
                    [valid j_nbr k_nbr j_spr k_spr] = get_indices(leaflet, j, k, j_nbr_tmp, k_nbr_tmp); 
                    
                    if valid 
                        X_nbr = X_current(:,j_nbr,k_nbr); 
                        
                        tension = tension_with_reference(X, X_nbr, R_v(j_spr,k_spr), k_v(j_spr,k_spr), leaflet); 
                        F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 
                    end 
                    
                end 

                % current node has a chordae connection
                if chordae_idx(j,k).tree_idx
                    
                    tree_idx = chordae_idx(j,k).tree_idx; 

                    [m N_chordae] = size(chordae(tree_idx).C);

                    % index in current free edge array 
                    i = chordae_idx(j,k).leaf_idx;
                    
                    % index that free edge would have if on tree
                    % remember that leaves are only in the leaflet
                    leaf_idx = chordae_idx(j,k).leaf_idx + N_chordae;

                    % then take the parent index of that number in chordae variables
                    idx_chordae = floor(leaf_idx/2);

                    X_nbr = chordae(tree_idx).C(:,idx_chordae);
                    
                    tension = tension_with_reference(X, X_nbr, chordae(tree_idx).R_free_edge(i), chordae(tree_idx).k_free_edge(i), leaflet);

                    F_tmp = F_tmp + tension * (X_nbr-X)/norm(X_nbr-X); 

                end 
                
                
                F_leaflet(:,j,k) = F_tmp;

            end
        end 
    end
    

    % chordae internal terms 
    for tree_idx = 1:num_trees
        
        C = chordae(tree_idx).C; 
        [m N_chordae] = size(C);         
        F_chordae(tree_idx).C = zeros(size(C));  
        
                
        for i=1:N_chordae

            left   = 2*i; 
            right  = 2*i + 1;
            parent = floor(i/2); 

            for nbr_idx = [left,right,parent]

                % get the neighbors coordinates, reference length and spring constants
                % routine handles unpacking and pulling correct constants 
                [nbr R_nbr k_val] = get_nbr_chordae(leaflet, i, nbr_idx, tree_idx); 
                
                tension = tension_with_reference(C(:,i), nbr, R_nbr, k_val, leaflet); 
                
                tension_by_tangent = tension * (nbr - C(:,i)) / norm(nbr - C(:,i));  

                F_chordae(tree_idx).C(:,i) = F_chordae(tree_idx).C(:,i) + tension_by_tangent;

            end 

        end 
    end 

    F = linearize_internal_points(leaflet, F_leaflet, F_chordae); 
    
end 



