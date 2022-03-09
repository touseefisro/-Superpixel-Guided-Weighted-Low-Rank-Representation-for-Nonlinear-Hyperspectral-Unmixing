function [SAD,AVG_SAD, MSE] = angle_gbm( W,A )

[M,P]=size(A);
for i = 1:P
      SAD_cor(i) = W(:, i)'*A(:, i)/(norm(W(:, i))*norm(A(:, i)));
end
SAD =acos(SAD_cor);

AVG_SAD= mean(SAD);
MSE=norm(W-A,'fro')/(sqrt(P*M));

end

