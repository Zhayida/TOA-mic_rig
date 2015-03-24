function nbrSols=testNbrSols_M2(eqs)

% addpath ../../../../../multipol/
% addpath ../../../../../polybase/matlab/


filename = ['testIdeal''.m2'];



%% take the eqs and replace all coeffecienst with numbers from Z/30097

for kk=1:length(eqs)
    mons=monomials(eqs(kk));
 %remaking the monomials to a mutip
    coeffsZMod=floor(30097*rand(size(mons,2),1));
    coeffsZMod=coeffsZMod';
    eqs(kk)=multipol(coeffsZMod,mons);
    
end


%% pick out all sub-matrices of 4 by 4


mm = eqs2m2(eqs,filename);

% run m2 script
eval(['!m2<',filename]);

end