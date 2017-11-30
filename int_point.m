clc

% TODO: Amat = randn(50,200) * 1/ sqrt(50);
%       xtrue 3-sparse vector


%% Amat is m x n where n > m
m = 3;
n = 8;
rho = 0.8;
mu = 10;
xtplus = zeros(n,1);
xtminus = zeros(n,1);
eps = 0.01;
c = ones(n*2,1);
e = ones(n*2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H2mat = [1 -1; 1 1];
H3mat = [H2mat -1*H2mat; H2mat H2mat];
H4mat = [H3mat -1*H3mat; H3mat H3mat];

Amat = H4mat([2,5,8],:);

Acon = [Amat -Amat];

xtrue = [0 0 3 0 0 0 0 -3]';
b = Amat * xtrue;

% xk = (Amat') * inv(Amat*Amat') * b;
[U,S,V] = svd(Amat);
xk = xtrue + 2*V(:,4) + 3*V(:,5); 
for i=1:n
    if xk(i) > 0
        xtplus(i) = xk(i);
    
    elseif xk(i) < 0
        xtminus(i) = -1 * xk(i);
    end
end


xk = [xtplus; xtminus];

alpha = 1;

for k=1:30
    if k>1
        mu = mu * rho;
    end
   for j=1:200
       
       X = diag(xk) + eps * eye(2*n,2*n);
       lambda = inv(Acon * (X * X) * Acon') * (Acon * (X * X) * c - mu * Acon * X * e);
       pb = X * ones(2*n,1) + (1 / mu) * X * X * (Acon' * lambda - c);
       
%        gradient = c - mu * inv(X) * e; 
%        hessian = mu * inv(X*X); 
       
%        ss = pb;
%        if gradient' * ss < 0
%            ss = (-1) * ss;
%        end
%        
%        numr = (gradient' * ss);
%        denom = (ss' * hessian * ss);
%        alpha = numr / denom;
%       
       xk = xk + alpha * pb;
   
   end
end




for j=1:2*n
    if abs(xk(j) < eps)
        xk(j) = 0;
    end
end
xk

% newx = [xk(1:n); -xk(n+1:2*n)]

% first = xk(1:n);
% second = xk(n+1:2*n);
% 
% newxk = first - second;

% 
% Acon * xk - b;


   
