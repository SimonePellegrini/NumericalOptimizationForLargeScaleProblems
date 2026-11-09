function grad = bandedTrigonometric_FDGrad(x,h,relative)
    n = length(x);
    if relative
        h = abs(x)*h;
    else
        h = ones(n,1)*h;
    end
    
    grad = ((1:n)'.*sin(x).*sin(h)+2.*cos(x).*sin(h))./h;
    grad(n) = (n*sin(x(n)).*sin(h(n))-(n-1)*cos(x(n)).*sin(h(n)))/h(n);
end

