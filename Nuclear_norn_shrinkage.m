function [ J] = Nuclear_norn_shrinkage( X, tau)
[U,sigma,V] = svd(X,'econ');
sigma = diag(sigma);
svp = length(find(sigma>tau));
if svp>=1
    sigma = sigma(1:svp)-tau;
else
    svp = 1;
    sigma = 0;
end
J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
