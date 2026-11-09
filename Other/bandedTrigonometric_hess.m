function H = bandedTrigonometric_hess(x)
    
    x = x(:);
    n = length(x);
    i_vec = (1:n)';
    

    diag_vals = i_vec .* cos(x) - 2 * sin(x);
    
    diag_vals(n) = n * cos(x(n)) + (n-1) * sin(x(n));
   
    H = spdiags(diag_vals, 0, n, n);
end