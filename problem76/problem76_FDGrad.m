function g = problem76_FDGrad(x,h,relative)
    n = length(x);
    if relative
        h = abs(x)*h;
    else
        h = ones(n,1)*h;
    end
    x_prev = [x(n);x(1:n-1)];
    x_next = [x(2:n);x(1)];
    g = (-10*x_prev.*x+x.^3+x.*h.^2+50*x-5.*x_next.^2)/50;
end

