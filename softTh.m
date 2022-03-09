%% Simple soft-thresholding
function X=softTh(B,lambda)
X=sign(B).*max(0,abs(B)-lambda);
end