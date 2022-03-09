function [SRE] = Evaluate_SRE( A, A1)

SRE=10*log10(norm(A,'fro')^2/norm(A-A1,'fro')^2);
