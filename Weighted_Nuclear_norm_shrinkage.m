% function [ J] = Weighted_Nuclear_norm_shrinkage( X, tau)
% [U,sigma,V] = svd(X,'econ');
% sigma = diag(sigma);
% svp = length(find(sigma>tau));
% if svp>=1
%     sigma = sigma(1:svp)-tau;
% else
%     svp = 1;
%     sigma = 0;
% end
% J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
% 
% 



function [X] = Weighted_Nuclear_norm_shrinkage(B,lambda)
%weighted singular value threshold function SVT
% The proximal operator of the nuclear norm of a matrix
% 
% min_X lambda*||X||_*+0.5*||X-B||_F^2
%
% version 1.0 - 06/09/2016
%
% Written by Touseef Ahmad; SAC-ISRO
% 

[U,W,V] = mySVD(B);
W = diag(W);
svp = length(find(W>lambda));
if svp>=1
    b = 1./(W(1:svp)+0.00001);
    W1 = max(W(1:svp)-lambda*b,0);
    X = U(:,1:svp)*diag(W1)*V(:,1:svp)';
    nuclearnorm = sum(W1);
else
    X = zeros(size(B));
    nuclearnorm = 0;
end

