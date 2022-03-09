%% ===============================Parameter analysis===========================================
clc;
clear all;
close all;
addpath(genpath(cd));

%% ============================Generate simulated data set============================
load USGS_pruned_10_deg.mat
M = 10;
E = B(:,randperm(62,M));
clearvars B;
[C, M] = size(E);
Row = 5;
Col = 5;
P = Row*Col;
Rank = 5;
Sparsity = 0.7;
A = Generate_abundance(M, P, Rank, Sparsity);
imagesc(A)
rank(A)

SNR = 30;
Nonlinear = 1;
Sparse = 1;
Stripes_num = 1;
[Y_orig, F, B, S, Gauss] = getSynData_type(E, A, Row, Col, C, SNR, Nonlinear, Sparse, Stripes_num);

Y = reshape(Y_orig, Row*Col, C)';


%% ============================RGBM-LRR==============================================
disp('RGBM-LRR');
mu = 1e-2;
tol = 1e-6;
maxiter = 100;
lambda = [0 10.^(-6:1:-1)];
alpha = [0 10.^(-6:1:-1)];

no_l1 = size(lambda,2);
no_l2 = size(alpha,2);

for j1=1:no_l1
    for j2 = 1:no_l2
        tic;
        Results5{j1,j2} = RGBM_LRR(Y, E, lambda(j1), alpha(j2), mu, tol, maxiter);
        t5(j1,j2) = toc;
        SRE5(j1,j2) = Evaluate_SRE(A, Results5{j1,j2}.A);
    end
end
[Results5_max, SRE5_max, Err5_max] = Get_max(Results5, SRE5, A);




%% ============================RGBM-SS-LRR==============================================
disp('RGBM-SS-LRR');
mu = 1e-2;
tol = 1e-6;
maxiter = 500;
K = 10; %Number of Super Pixels
lambda = [0 10.^(-6:1:-1)];
alpha = [0 10.^(-6:1:-1)];

no_l1 = size(lambda,2);
no_l2 = size(alpha,2);

Y_3D = reshape(Y', Row, Col, C);
for j1=1:no_l1
    %j1
    for j2 = 1:no_l2
        tic
        [Results2{j1,j2}] = RGBM_SS_LRR(Y_3D, E, K, lambda(j1), alpha(j2), mu, tol, maxiter);
        t2(j1,j2) = toc;
        [SAD2(j1,j2,:),AVG_SAD2(j1,j2), RE2(j1,j2)] = angle_gbm( Y, Results2{j1,j2}.Y);
        
          SRE2(j1,j2) = Evaluate_SRE(A, Results2{j1,j2}.A);
        
    end
end
[Results2_max, SRE2_max, Err2_max] = Get_max(Results2, SRE2, A);




%%%%%%%%%%%%%%%.... WEIGHTED LRR & SUPER-PIXELS GUIDED WEIGHTED LRR
%% ============================RGBM-WLRR==============================================
disp('RGBM-WLRR');
mu = 1e-2;
tol = 1e-6;
maxiter = 500;
lambda = [0 10.^(-6:1:-1)];
alpha = [0 10.^(-6:1:-1)];

no_l1 = size(lambda,2);
no_l2 = size(alpha,2);

for j1=1:no_l1
    for j2 = 1:no_l2
        tic;
        Results5_w{j1,j2} = RGBM_WLRR(Y, E, lambda(j1), alpha(j2), mu, tol, maxiter);
        t5_w(j1,j2) = toc;
        
       [SAD5_w(j1,j2,:),AVG_SAD5_w(j1,j2), RE5_w(j1,j2)] = angle_gbm( Y, Results5_w{j1,j2}.Y);

        SRE5_w(j1,j2) = Evaluate_SRE(A, Results5_w{j1,j2}.A);
    end
end
[Results5_w_max, SRE5_w_max, Err5_w_max] = Get_max(Results5_w, SRE5_w, A);



%% ============================RGBM-SG-WLRR==============================================
disp('RGBM-SS-WLRR');
mu = 1e-2;
tol = 1e-6;
maxiter = 500;
K_s = 10; %Number of Super Pixels
lambda = [0 10.^(-6:1:-1)];
alpha = [0 10.^(-6:1:-1)];

no_l1 = size(lambda,2);
no_l2 = size(alpha,2);

% Y_3D = reshape(Y', Row, Col, C);
for j1=1:no_l1
    %j1
    for j2 = 1:no_l2
        tic
        [Results2_w{j1,j2}] = RGBM_SG_WLRR(Y_3D, E, K_s, lambda(j1), alpha(j2), mu, tol, maxiter);
        t2_w(j1,j2) = toc;
        [SAD2_w(j1,j2,:),AVG_SAD2_w(j1,j2), RE2_w(j1,j2)] = angle_gbm( Y, Results2_w{j1,j2}.Y);
        
        SRE2_w(j1,j2) = Evaluate_SRE(A, Results2_w{j1,j2}.A);
        
    end
end
[Results2_w_max, SRE2_w_max, Err2_w_max] = Get_max(Results2_w, SRE2_w, A);









