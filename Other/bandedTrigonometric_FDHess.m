function hess = bandedTrigonometric_FDHess(x,h)
    n=length(x);
    hess = (1-cos(h)).*(2*(1:n)'.*cos(x)-4.*sin(x))/h^2; 
    hess(n) = (1-cos(h)).*(2*n.*cos(x(n))+2*(n-1)*sin(x(n)))/h^2; 
    hess = spdiags(hess,0,n,n);
end

