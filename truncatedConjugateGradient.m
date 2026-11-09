function [z_j,num_iter] = truncatedConjugateGradient(gradfk,hessian,max_iter)
   
    % INPUT:
    %   gradfk   - The gradient of the objective function at the current iterate x_k.
    %   hessian  - The exact Hessian matrix of the objective function evaluated at x_k.
    %   max_iter - Maximum number of inner Conjugate Gradient iterations allowed.
    
    % OUTPUT:
    %   z_j      -  The computed search direction 
    %   num_iter -  The actual number of inner iterations performed.

        z_j = zeros(length(gradfk),1);
        c_k = -gradfk;
        r_j = -c_k;
        d_j = -r_j;
        norm_c_k = norm(c_k);
        
        %eta_k selected to guarantee a superlinear/quadratic
        %rate of convergence
        eta_k = min(0.5,sqrt(norm_c_k));
        num_iter = 0;   
        
        %check the stopping criterion
        while norm(r_j)>eta_k*norm_c_k && num_iter<max_iter
            %pre compute some products to avoid computing
            %some things multiple times
            B_kd_j  = hessian*d_j;
            d_jTB_kd_j = d_j'*B_kd_j;
            if d_jTB_kd_j>0
                %we trust the matrix is positive definite
                r_jTr_j = r_j'*r_j;
                %compute the alpha that minimizes J(x)
                alpha_j = (r_jTr_j)/(d_jTB_kd_j);

                %compute the new point
                z_j = z_j + alpha_j*d_j;
                
                %compute the new direction
                r_jp1 = r_j + alpha_j*B_kd_j;
                beta_k = (r_jp1'*r_jp1)/(r_jTr_j);
                d_j = -r_jp1+beta_k*d_j;

                %update the variables
                r_j = r_jp1;
                num_iter = num_iter+1;
            elseif num_iter==0
                %set the computed direction equal to the negative gradient
                num_iter = num_iter+1;
                z_j = -gradfk;
                break;
            else
                num_iter = num_iter+1;
                break;
            end
        end
end

