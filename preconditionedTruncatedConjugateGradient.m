function [z_j,num_iter] = preconditionedTruncatedConjugateGradient(gradfk,hessian,L,max_iter)
        z_j = zeros(length(gradfk),1);
        c_k = -gradfk;   
        r_j = gradfk;
        w_j = L\r_j;
        y_j = L'\w_j;   
        d_j = -y_j;
        %eta_k selected to guarantee a superlinear/quadratic
        %rate of convergence
        norm_c_k = norm(c_k);
        eta_k = min(0.5,sqrt(norm_c_k));
        num_iter = 0;   
        
        %check the stopping criterion
        while norm(r_j)>eta_k*norm_c_k && num_iter<max_iter
            %pre compute some products to avoid computing
            %some things multiple times
            Hd_j  = hessian * d_j;
            d_jTHd_j = d_j'*Hd_j;
            if d_jTHd_j>0
                %we trust the matrix is positive definite
                r_jTy_j = r_j'*y_j;
                %compute the alpha that minimizes J(x)
                alpha_j = (r_jTy_j)/(d_jTHd_j);

                %compute the new point
                z_j = z_j + alpha_j*d_j;
                
                %compute the new direction
                r_jp1 = r_j + alpha_j*Hd_j;
                w_jp1 = L\r_jp1;
                y_jp1 = L'\w_jp1;

                beta_j = (r_jp1'*y_jp1)/(r_jTy_j);
                d_j = -y_jp1+beta_j*d_j;

                
                %update the variables
                r_j = r_jp1;
                y_j = y_jp1;
                num_iter = num_iter+1;
            elseif num_iter==0
                num_iter = num_iter+1;
                z_j = c_k;
                return
            else
                num_iter = num_iter+1;
                return
            end
        end       
end

