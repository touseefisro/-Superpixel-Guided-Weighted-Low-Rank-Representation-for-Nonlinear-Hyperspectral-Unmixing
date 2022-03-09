function A = Generate_abundance(M, P, Rank, Sparsity)
if Rank < 0 || Rank > P || Sparsity < 0 || Sparsity >= 1
    error('The Rank or Sparsity is wrong');
end
A = rand(M,P);
Res = mod(P, Rank);
Num = floor(P/Rank);
if Res ~= 0
    for i=1:Rank-1
        I = randperm(M,round(Sparsity*M));
        A(:,1+(i-1)*Num:i*Num) = rand(M,1)*rand(1,Num);
        A(I,1+(i-1)*Num:i*Num) = 0;
    end
    I = randperm(M,round(Sparsity*M));
    A(:,1+(Rank-1)*Num:end) = rand(M,1)*rand(1,Num + Res);
    A(I,1+(Rank-1)*Num:end) = 0;
else
    for i=1:Rank
        I = randperm(M,round(Sparsity*M));
        A(:,1+(i-1)*Num:i*Num) = rand(M,1)*rand(1,Num);
        A(I,1+(i-1)*Num:i*Num) = 0;
    end
end
