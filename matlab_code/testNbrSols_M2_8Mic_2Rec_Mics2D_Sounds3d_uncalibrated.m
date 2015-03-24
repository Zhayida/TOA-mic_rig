% determinant constraints

addpath ../../../../../multipol/
addpath ../../../../../polybase/matlab/
clear
%%
% m = 9;n = 5; % (0,1)  % 70  eqs % 5 unknowns
% m = 8;n = 5; % (1,5)  
% m = 7;n = 5; % (2,10) 
% m = 7;n = 6; % (0,5)  % 75  eqs % 6 unknowns
% m = 6;n = 8; % (0,14) % 175 eqs 
% m = 6;n = 6; % (2,30)
% m = 6;n = 7; % (1,35)

filename = ['tester_',num2str(m),'_',num2str(n),'.m2'];

if ~exist(filename,'file') && 1
    
    
    %% rank-3 constraints on sub-determinants
    D = round(10*abs(randn(m,n)));
    
    V = eye(n);
    for i = 1:n;
        tv(i) = multipol(1,V(:,i));
    end
    one = multipol(1,zeros(n,1));
    zero =multipol(0,zeros(n,1));
    
    Dv = D*one;
    clear Tv;
    for i = 2:m
        for j = 2:n
            Tv(i-1,j-1) = Dv(i,j)^2 - 2*Dv(i,j)*tv(j) - ...
                Dv(1,j)^2 + 2*Dv(1,j)*tv(j)     - ...
                Dv(i,1)^2 + 2*Dv(i,1)*tv(1)     + ...
                Dv(1,1)^2 - 2*Dv(1,1)*tv(1) ;
        end
    end
        
    
    %% pick out all sub-matrices of 4 by 4
    rp = nchoosek([1:m-1],4);
    cp = nchoosek([1:n-1],4);
    
    clear eqs;
    k = 1;
    for i = 1:size(rp,1)
        for j = 1:size(cp,1)
            CC = Tv(rp(i,:),cp(j,:));
            eqs(k) = detv(CC);
            k = k+1;
        end
    end
    
    mm = eqs2m2(10*eqs,filename);
end
% run m2 script
eval(['!m2<',filename]);