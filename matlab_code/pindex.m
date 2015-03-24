function [P,pvec] = pindex(n);
% generate pairwise indices in a matrix
% e.g. n = 3, 
% P = [ 1 2 3;
%       2 4 5;
%       3 5 6];

k = 0;
P = zeros(n);
for i = 1:n
    for j = i:n
        k = k+1;
        P(i,j) =k;
        pvec(k,:) = [i,j];
        P(j,i) = k;
    end
end
