function hess = bandedTrigonometric_hessGivenGrad(x,h,relative)
n = length(x);
    if relative
        h = h*abs(x);
    else
        h = ones(n,1)*h;
    end
    diag = (-2*sin(x).*sin(h)+(1:n)'.*cos(x).*sin(h))./h;
    diag(n)= ((n-1)*sin(x(n))*sin(h(n))+n*cos(x(n))*sin(h(n)))/h(n);
    hess = spdiags(diag,0,n,n);
end

