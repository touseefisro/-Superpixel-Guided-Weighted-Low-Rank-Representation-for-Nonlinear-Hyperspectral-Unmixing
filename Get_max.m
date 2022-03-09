function [Results_max, SRE_max, Err_max, Index] = Get_max(Results, SRE, A)
[~,ind ] = max(SRE(:));
[Row, Col] = ind2sub(size(SRE),ind);
Results_max = Results{Row,Col};
SRE_max = SRE(Row,Col);
Err_max = abs(A-Results_max.A);
Index = [Row, Col];

