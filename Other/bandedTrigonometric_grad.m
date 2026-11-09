function g = bandedTrigonometric_grad(x)
    
    n = length(x);
    i_vec = (1:n)';
    
    g = i_vec .* sin(x) + 2 * cos(x);

    g(n) = n * sin(x(n)) - (n-1) * cos(x(n));
end