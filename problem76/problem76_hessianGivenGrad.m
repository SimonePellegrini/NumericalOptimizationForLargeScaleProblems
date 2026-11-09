function H = problem76_hessianGivenGrad(x,h,relative)


n = length(x);
d = zeros(3,n);


if relative
    h = abs(x)*h;
else
    h = ones(length(x),1)*h;
end

for c = (1:3)
    e = zeros(n,1);
    idx = c+1 : 3 : n-1;
    e(idx) = 1;
    if relative
        h_eff = zeros(n,1);

        h_eff(c:3:n-1) = h(c+1:3:n);

        for i = 2:n
            if h_eff(i) == 0
                h_eff(i) = h_eff(i-1);
            end
        end

        d(c,:) = (problem76_grad(x+h.*e)-problem76_grad(x-h.*e))./(2*h_eff);
    else
        d(c,:) = (problem76_grad(x+e.*h)-problem76_grad(x-e.*h))./(2*h);
    end
end

e = zeros(n,1);
e(n) = 1;
last_col = (problem76_grad(x+e.*h)-problem76_grad(x-e.*h))/(2*h(n));

e(n)=0;
e(1)=1;
first_col = (problem76_grad(x+e.*h)-problem76_grad(x-e.*h))/(2*h(1));

mainDiag = zeros(n,1);
upperDiag = zeros(n,1);
lowerDiag = zeros(n,1);

mainDiag(1) = first_col(1);
mainDiag(2:3:n) = d(1,2:3:n);
mainDiag(3:3:n) = d(2,3:3:n);
mainDiag(4:3:n) = d(3,4:3:n);
mainDiag(n) = last_col(n);


upperDiag(2:3:n) = d(1,1:3:n-1);
upperDiag(3:3:n) = d(2,2:3:n-1);
upperDiag(4:3:n) = d(3,3:3:n-1);
upperDiag(n) = last_col(n-1);

lowerDiag(1) = first_col(2);
lowerDiag(2:3:n-1) = d(1,3:3:n);
lowerDiag(3:3:n-1) = d(2,4:3:n);
lowerDiag(4:3:n-1) = d(3,5:3:n);
lowerDiag(n) = last_col(n);


H = spdiags([lowerDiag mainDiag upperDiag], [-1 0 1], n, n);


H(1,n) = last_col(1);
H(n,1) = first_col(n);


end


