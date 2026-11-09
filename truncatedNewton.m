function [gradfk_norm,fk,errorSeq,x_k,k,xseq] = ...
    truncatedNewton(x_0, f, gradf,hessianf, ...
    max_iter, tolgrad, c1, rho, btmax,pcg_max_iter, precondition)
    
    % INPUT:
%   x_0      - Starting point of the iterative method.
%   f        - Function handle: Objective function to be minimized.
%   gradf    - Function handle: Exact gradient of the objective function.
%   hessianf - Function handle: Exact Hessian matrix of the objective function.
%   max_iter - Maximum number of outer iterations allowed.
%   tolgrad  - Tolerance for the gradient norm stopping criterion.
%   c1       - Armijo condition parameter.
%   rho      - Step size reduction factor (0 < rho < 1) for backtracking.
%   btmax    - Maximum number of backtracking steps per iteration.
  
% OUTPUT:
%   gradfk_norm - Scalar: Euclidean norm of the gradient at the final iteration.
%   fk          - Scalar: Value of the objective function at the computed minimum.
%   errorSeq    - Vector: Sequence of step lengths (norm(x_{k+1} - x_k)) for each iteration.
%   x_k         - Vector: The final approximation of the minimizer found by the algorithm.
%   k           - Integer: The total number of outer iterations performed.
%   xseq        - Matrix: History of the iterations, where the i-th row is the point x at iteration i.

    %computes the value used in the armijo condition
    farmijo = @(fk, alpha, c1_gradfk_pk) fk + alpha * c1_gradfk_pk;
    
    %set up the initial variables
    x_k = x_0;
    fk = f(x_k);
    xseq = x_k(:)';
    errorSeq = zeros(max_iter,1);
    total_inner_iterations = 0;
    %computes the gradient to build the second order Taylor's expansion
    gradfk = gradf(x_k);
    %computes the hessian to build the second order Taylor's expansion
    hessian = hessianf(x_k);

    %computes the norm of the gradient
    gradfk_norm = norm(gradfk); 
    
    %variable that counts the number of iterations
    k = 0;
    %main loop
    while k < max_iter && gradfk_norm >= tolgrad
        if precondition
            try 
                L = ichol(hessian);
                [z_j,inn_iterations] = preconditionedTruncatedConjugateGradient(gradfk,hessian,L,pcg_max_iter);
            catch
                [z_j,inn_iterations] = truncatedConjugateGradient(gradfk,hessian,pcg_max_iter); 
            end
        else
            [z_j,inn_iterations] = truncatedConjugateGradient(gradfk,hessian,pcg_max_iter);
        end
        total_inner_iterations = total_inner_iterations + inn_iterations;
        %set z as the new direction for the inexact linesearch
        pk = z_j;
        k = k + 1;
        alpha = 1;

        
         %computes phi'(0)
        c1_gradfk_pk = c1 * gradfk' * pk;
        bt = 0;

        %Update xnew and fnew with the new point computed
        xnew = x_k+alpha*pk;
        fnew = f(xnew);
        
        %Armijo + BackTracking to compute the direction of the search
        while bt < btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
            % Reduce the value of alpha
            alpha = rho * alpha;
            % Update xnew and fnew w.r.t. the reduced alpha
            xnew = x_k + alpha * pk;
            fnew = f(xnew);
            % Increase the counter by one
            bt = bt + 1;
        end
        

        %if we didn't find a new point satisfying the armijo condition then
        %end the loop
        if bt == btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
            display("Armijo condition not satisfied!")
            break
        end

        errorSeq(k) = norm(x_k-xnew);
        x_k = xnew;

        % update the value of the function, the gradient, the hessian and
        % the norm of the gradient.
        xseq(end+1,:) = x_k(:)';
        fk = f(x_k);
        gradfk = gradf(x_k);
        %computes the hessian to build the second order Taylor's expansion
        hessian = hessianf(x_k);
        
        gradfk_norm = norm(gradfk);
    end
    errorSeq = errorSeq(1:k);
    
end

 