function Results = RGBM_LRR(Y, E, lambda, alpha, mu, tol, maxiter)
%% Check input
if nargin < 2
   error('The input data should be at least two!');
end
[C, P] = size(Y);
M = size(E, 2);
%% Initializing optimization variables
A = Fcls(Y, E);
F = Bilinear_endmember( E );
B = zeros(size(F,2), size(Y, 2));
S = zeros(size(Y));

%Initialization of auxiliary matrices V1,V2,V3,V4
V1 = A;
V2 = A;
V3 = B;

%Initialization of Lagranfe Multipliers
D1 = V1*0;
D2 = V2*0;
D3 = V3*0;

%primal residual 
res_p = inf;

%dual residual
res_d = inf;

%indication variable
mu_changed = 0;

%error tolerance
tol1 = sqrt((3*M + C)*P)*tol;

%current iteration number
k=1;

%% Main iteration
while (k < maxiter) && ((abs (res_p) > tol1) || (abs (res_d) > tol1))
    if mod(k, 10) == 1
        %max(res_p,res_d)
        V10 = V1;
        V20 = V2;
        V30 = V3;
    end
    V1 = Nuclear_norn_shrinkage(A + D1, lambda/mu);
    A = inv(E'*E + 2*mu*eye(size(E,2)))*(E'*(Y - F*B - S) + mu*(V1 - D1 + V2 - D2));
    V2 = max(A + D2, 0);
    B = inv(F'*F + mu*eye(size(F,2)))*(F'*(Y - E*A - S)+ mu*(V3 - D3));
    V3 = min( max(B + D3, 0), Bilinear_abubdance( A ));
    S = softTh(Y - E*A - F*B, alpha);
    
    % Lagrange multipliers update
    D1 = D1 - (V1 - A);
    D2 = D2 - (V2 - A);
    D3 = D3 - (V3 - B);
    
    % update mu so to keep primal and dual residuals whithin a factor of 10
    if mod(k, 10) == 1
        %primal residual
        res_p = norm([V1; V2; V3;] - [A; A; B;], 'fro');
        %dual residual
        res_d = mu*norm([V1; V2; V3;] - [V10; V20; V30;], 'fro');
        % update mu
        if res_p > 10*res_d
            mu = mu*2;
            D1 = D1/2;
            D2 = D2/2;
            D3 = D3/2;
            mu_changed = 1;
        elseif res_d > 10*res_p
            mu = mu/2;
            D1 = D1*2;
            D2 = D2*2;
            D3 = D3*2;
            mu_changed = 1;
        end
        if  mu_changed
            mu_changed = 0;
        end 
    end
    k = k + 1;
end
% if k == maxiter
%     display('Maximum iteration has been reached!');
% end
A(A<1e-4) = 0;
B(B<1e-4) = 0;
Results.A = A;
Results.F = F;
Results.B = B;
Results.S = S;
Results.Y=E*A+F*B;
end

function F = Bilinear_endmember( E )
M = size(E,2);
num1 = 1;
for k=1:M-1
    for j=k+1:M
        F(:,num1) = E(:,k).*E(:,j);
        num1=num1+1;
    end
end
end

function B = Bilinear_abubdance( A )
M = size(A,1);
num1 = 1;
for k=1:M-1
    for j=k+1:M
        B(num1,:) = A(k,:).*A(j,:);
        num1=num1+1;
    end
end
end