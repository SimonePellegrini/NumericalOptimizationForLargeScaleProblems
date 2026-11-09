function G = problem76_grad(x)
    n = length(x);
    x_next = [x(2:n);x(1)];
    f = x - (x_next.^2) / 10;
    f_next = [f(n);f(1:n-1)];
    G = f-f_next.*x/5;
end