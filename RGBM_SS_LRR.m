function Results  = RGBM_SS_LRR( Y, E, K, lambda, alpha, mu, tol, maxiter)
%%  RGBM_LRR: Robust GBM Low-Rank-Represenation & Superpixel-segmentation

%% Obtain principle components
[M, N, C] = size(Y);
p = 1; % number of principal components
Y_2d=reshape(Y, M*N, C);
[X_pca] = pca(Y_2d');%, p);

X_pca=X_pca(:,1);
img = im2uint8(mat2gray(reshape(X_pca, M, N, p)));
%% Superpixel segmentation
labels = mex_ers(double(img), K);
% grey_img = im2uint8(mat2gray(Y(:,:,1)));
% [height width] = size(grey_img);
% [bmap] = seg2bmap(labels,width,height);
% bmapOnImg = img;
% idx = find(bmap>0);
% timg = grey_img;
% timg(idx) = 255;
% bmapOnImg(:,:,2) = timg;
% bmapOnImg(:,:,1) = grey_img;
% bmapOnImg(:,:,3) = grey_img;
% 
% figure;
% imshow(bmapOnImg,[]);
Results_segment= seg_im_class(Y, labels);
Num=size(Results_segment.Y,2);
%% LRR
for i=1:Num
    Results_temp= RGBM_LRR(Results_segment.Y{1,i}', E, lambda, alpha, mu, tol, maxiter);
    A(:, Results_segment.index{1,i}) = Results_temp.A;
    B(:, Results_segment.index{1,i}) = Results_temp.B;
    S(:, Results_segment.index{1,i}) = Results_temp.S;
end
Results.A = A;
Results.F = Results_temp.F;
Results.B = B;
Results.S = S;
Results.Y=E*A+Results_temp.F*B;
end
