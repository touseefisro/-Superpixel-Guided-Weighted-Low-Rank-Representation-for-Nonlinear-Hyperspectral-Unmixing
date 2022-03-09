function [Y, F, B, S, Gauss] = getSynData_type(E, A, Row, Col, C, SNR, Nonlinear, Sparse, Stripes_num)

X = reshape((E*A)', Row, Col, C);
[M, N, C] = size(X);
%add Gauss noise
if SNR == 0
    Gauss = 0;
else
    Gauss = addnoise_hyperspectral_gauss(X, SNR);
end
% Nonlinear
[F, B] = GBM(E, A, Nonlinear);
% add Sparse noise
if Sparse==0 
    S = 0;
else
    S = addnoise_hyperspectral_sparse(X + reshape((F*B)', Row, Col, C) + Gauss, Stripes_num);
end


Y = X + reshape((F*B)', Row, Col, C) + S + Gauss;
%Y(Y<0) = 0;
end


function [F, B] = GBM(E, A, type)
    M = size(E,2);
    Gamma=unifrnd(0,1,M,M);
    num1 = 1;
    for k=1:M-1
        for j=k+1:M
            F(:,num1) = E(:,k).*E(:,j);
            if type == 0
                B(num1,:) = 0.*A(k,:).*A(j,:);
            else
                B(num1,:) = Gamma(k,j).*A(k,:).*A(j,:);
            end
            num1=num1+1;
        end
    end
end

function [ Gauss ] = addnoise_hyperspectral_gauss(X, SNR)
[M,N,B] = size(X);
for i=1:B
    Y(:,:,i) = awgn(X(:,:,i), SNR,'measured');
end
Gauss = Y - X;
end

function [ S ] = addnoise_hyperspectral_sparse(X, Stripes_num)
[M,N,B] = size(X);
Y = X;
Impulse_range = [20:30];
Stripes_range = [80:90];
for j=Impulse_range(1):Impulse_range(end)
    Y(:,:,j) = imnoise( X(:,:,j),'salt & pepper',0.1);
end


for k=Stripes_range(1):Stripes_range(end)
    temp = randperm(N);
    Select = temp(1:Stripes_num);
    for l = 1:Stripes_num
        Y(:,Select(l),k) = 0;
    end
end
S = Y - X;
end
