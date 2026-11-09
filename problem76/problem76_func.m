function F = problem76_func(x)

    n = length(x);
    x_next = [x(2:n); x(1)];
    f = x - (x_next.^2) / 10;
    F = 0.5 * sum(f.^2);
    
end