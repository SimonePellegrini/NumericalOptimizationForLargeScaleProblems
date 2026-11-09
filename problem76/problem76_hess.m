function H = problem76_hess(x)

n = length(x);

x_next = [x(2:n); x(1)];
f = x - (x_next.^2) / 10;
f_prev = [f(n); f(1:n-1)];

mainDiag  = 1 + (x.^2)/25 - f_prev/5;
outerDiag = -x/5;
H = spdiags([[outerDiag(2:end);0], mainDiag, [0;outerDiag(2:end)]], [-1, 0, 1], n, n);

if n>2
    H(1,n) = -x(1)/5;
    H(n,1) = -x(1)/5;
else
    H(1,n) = -(x(n)+x(1))/5;
    H(n,1) = -(x(n)+x(1))/5;
end

end
