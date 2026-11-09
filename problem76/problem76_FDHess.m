function hess = problem76_FDHess(x,h)
    n = length(x);
    x_prev = [x(n);x(1:n-1)];
    diag = (-20*x_prev+6*x.^2+h^2+100)/100;
    lower_diag = (-2*x-h)/10;
    hess = spdiags([[lower_diag(2:end);0], diag, [0;lower_diag(2:end)]], [-1, 0, 1], n, n);
    if n~= 2
        hess(1,n) = lower_diag(1);
        hess(n,1)= lower_diag(1);
    else
        hess(1,n) = 1/5*(-x(1)-x(2)-h);
        hess(n,1)= 1/5*(-x(1)-x(2)-h);
    end
end

