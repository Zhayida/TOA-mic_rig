function [sols]=toa_rig_8_2(d)

%8 mics in 2d, 2 sounds in 3d

if ispc
addpath C:\Users\simonb\Documents\multipol
addpath .
end

%% Constants

nMics=8;




d2 = d.^2;

d11 = d2(1,1);
d1j = d2(1,:)-d2(1,1);
di1 = d2(:,1)-d2(1,1);
dij = d2-repmat(d2(:,1),1,2)-repmat(d2(1,:),nMics,1)+d2(1,1);
Rtilde = [0 dij(2:end,2)'];
Stilde = [0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% so far so good

%%% d11

%%%%%%%%%%%%%%%%
%%% d1j

% v1 = R(1,1)^2;
% v2 = R(1,2)^2;
% v3 = R(1,3)^2;
% v4 = R(1,4)^2;
% v5 = R(1,5)^2;
% v6 = R(1,6)^2;
% v7 = R(1,7)^2; 
% v8 = R(1,8)^2;
% u12 = R(1,1)*R(1,2);
% u34 = R(1,3)*R(1,4);
% u56 = R(1,5)*R(1,6);
% u78 = R(1,7)*R(1,8);



%%%%%%%%%%%%%%%%
%%% di1
%testa di(7)-
%di(8)-

B=zeros(10,1); %the right hand vector in our system of linear equations.
A=zeros(10,14); %10 eqs. (7 brom block B, and 3 from rigDist equations). 14 unknowns (10 originally, but then we parametrize u12, u34, u56, u78 to get linear equations)
B(1) = di1(2);
B(2) = di1(3);
B(3) = di1(4);
B(4) = di1(5);
B(5) = di1(6); %Unknowns are in order of ( h  b v1 v2 v3 v4 v5 v6 v7
B(6) = di1(7);
B(7) = di1(8);
A(1,1:14)=[Rtilde(2)^2 (-2*Rtilde(2)/(-2)) -1 1 0 0 0 0 0 0 0 0 0 0];
A(2,1:14)=[Rtilde(3)^2 (-2*Rtilde(3)/(-2)) -1 0 1 0 0 0 0 0 0 0 0 0];
A(3,1:14)=[Rtilde(4)^2 (-2*Rtilde(4)/(-2)) -1 0 0 1 0 0 0 0 0 0 0 0];
A(4,1:14)=[Rtilde(5)^2 (-2*Rtilde(5)/(-2)) -1 0 0 0 1 0 0 0 0 0 0 0];
A(5,1:14)=[Rtilde(6)^2 (-2*Rtilde(6)/(-2)) -1 0 0 0 0 1 0 0 0 0 0 0];
A(6,1:14)=[Rtilde(7)^2 (-2*Rtilde(7)/(-2)) -1 0 0 0 0 0 1 0 0 0 0 0];
A(7,1:14)=[Rtilde(8)^2 (-2*Rtilde(8)/(-2)) -1 0 0 0 0 0 0 1 0 0 0 0];
%xxx = [h b v1 v2 v3 v4 v5 v6 v7 v8 u12 u34 u56 u78]'; %yeah, derp.


%%%%%%%%%%%%%%%%
%% ridDist^2 -equations
% We now use the three equations ||r_1-r_2||=||r3-r4||=||r5-r6||=||r7-r8||


% ||r_1-r_2||-||r3-r4||=0 foloows here
 %norm(R(:,2)-R(:,1))
%U12^2
 %groud truth test, followed by Parametrization ground test. Each term
 %above and below should be equal.
% R(:,2)'*R(:,2)      -2*R(:,2)'*R(:,1) + R(:,1)'*R(:,1) - R(:,3)'*R(:,3)      + 2*R(:,3)'*R(:,4)                 - R(:,4)'*R(:,4) %ok
% v2 + Rtilde(2)^2*h  -2* u12           + v1             -  v3 - Rtilde(3)^2*h + 2*Rtilde(3)*Rtilde(4)*h + 2*u34  - v4-Rtilde(4)^2*h  %ok

B(8,1) = 0; %we take the two rig distances minus each other 
A(8,1:14)=[(Rtilde(2)^2-Rtilde(3)^2 + 2*Rtilde(3)*Rtilde(4) - Rtilde(4)^2) 0 1 1 -1 -1 0 0 0 0 -2 2 0 0];

%%%%%%%%%%%%%%%%
%% 

% ||r_1-r_2||-||r5-r6||=0 foloows here
% norm(R(:,6)-R(:,5)); %should be =rigDist

 %groud truth test, followed by Parametrization ground test. Each term
 %above and below should be equal.
% R(:,2)'*R(:,2)      -2*R(:,2)'*R(:,1)   + R(:,1)'*R(:,1)  - R(:,5)'*R(:,5)       + 2*R(:,5)'*R(:,6)                  - R(:,6)'*R(:,6) %ok
% 
% v2 + Rtilde(2)^2*h  -2* u12             + v1              - v5 - Rtilde(5)^2*h   + 2*Rtilde(5)*Rtilde(6)*h + 2*u56   - v6-Rtilde(6)^2*h  %ok

B(9) = 0;
A(9,1:14)=[(Rtilde(2)^2-Rtilde(5)^2 + 2*Rtilde(5)*Rtilde(6) - Rtilde(6)^2) 0 1 1 0 0 -1 -1  0 0 -2 0 2 0];


%%%%%%%%%%%%%%%%
%% 

% ||r_1-r_2||-||r7-r8||=0 foloows here
%  norm(R(:,8)-R(:,7)); %should be =rigDist

 %groud truth test, followed by Parametrization ground test. Each term
 %above and below should be equal.
% R(:,2)'*R(:,2)      -2*R(:,2)'*R(:,1)   + R(:,1)'*R(:,1)  - R(:,7)'*R(:,7)       + 2*R(:,7)'*R(:,8)                  - R(:,8)'*R(:,8) %ok

% v2 + Rtilde(2)^2*h  -2* u12             + v1              - v7 - Rtilde(7)^2*h   + 2*Rtilde(7)*Rtilde(8)*h + 2*u78   - v8-Rtilde(8)^2*h  

B(10) = 0;
A(10,1:14)=[(Rtilde(2)^2-Rtilde(7)^2 + 2*Rtilde(7)*Rtilde(8) - Rtilde(8)^2) 0 1 1 0 0 0 0 -1 -1  -2 0 0 2];


%% %%%%%%%%%%%%%%

%SiBu: chosen scientifcally using chaos theory, for the betterment of all
%mankind. ok, so: need to select ten cols that so that that part of A has
%full rank and low condition number. THis is the most stable choice i found
%after 50000 randozations over several constelaltions of ground truth.
extractedColsFromAInd= [4     5     6     7     8     9    10    12    13    14];

A1 = A(:,extractedColsFromAInd);
unknownsInd=setdiff(1:14,extractedColsFromAInd);
A2 = A(:,unknownsInd);
         
C1 = -A1\A2;
C2 = A1\B;
monvek(1,1)=multipol(1,[1;0;0; 0]);  %This is h
monvek(2,1)=multipol(1,[0;1;0; 0]);  %b
monvek(3,1)=multipol(1,[0;0;1; 0]);  %v1
monvek(4,1)=multipol(1,[0;0;0; 1]);  %u12


vu = C1*monvek+C2; %the other monomials expresad as al incomb of h, b,v u12
%Vu thus consists of the unknowns [v2 v3 v4 v5 v6 v7 v8 u34 u56 u78] inthat
%order, express in the unknowns h b v1 u12
%% THe equations to be solved by polysolve
% We here use the extra equations v1*v2=u12^2, v3*v4=u34^2, v5*v6=u56^2, v7*v8=u78^2
eqs(1,1)=monvek(3)*vu(1)-monvek(4)^2;  %v1*v2=u12^2
eqs(2,1)=vu(2)*vu(3)-vu(8)^2; %v3*v4=u34^2
eqs(3,1)=vu(4)*vu(5)-vu(9)^2; %v5*v6=u56^2
eqs(4,1)=vu(6)*vu(7)-vu(10)^2; %v7*v8=u78^2
%gt = xxx(unknownsInd);
%norm([evaluate(vu,gt) - xxx(extractedColsFromAInd)])/norm(xxx(extractedColsFromAInd))
%evaluate(eqs,gt)

%eqs2m2(100000*eqs,'eqs_8_2.m2')
%eval(['!m2<','eqs_8_2.m2'])
eqs2 = expandtodegree(eqs,[2 2 2 2],4);


% for i = 1:length(eqs)
%     [xx,yy] = polynomials2matrix(eqs(i));
%     yy = monvec2matrix(yy);
%     zz{i} = yy;
% end


% eqs2 = generate_equations(eqs,3);

settings.dim = 6;
%settings.basis_size = 14;
settings.basis_size = 10;
settings.action_variable = 1;
[psols stats] = polysolve(eqs2, settings);

sols=struct;
sols.R=cell(0);
sols.S=cell(0);


for kk = 1:size(psols,2);
    % For each solution try to backtrack and calculate all
    % unknown parameters
    onepsol = psols(:,kk);
    vuEval=evaluate(vu,onepsol);
    xxxsol = [onepsol(1:3) ; vuEval(1:7) ; onepsol(4) ; vuEval(8:10) ];  %recreating the unkowns [h b v1 v2 v3 v4 v5 v6 v7 v8 u12 u34 u56 u78]
    % xxx = [h b v1 v2 v3 v4 v5 v6 u12 u34 u56]';
    hh = xxxsol(1);
    ll = 1/sqrt(hh);
    bb = xxxsol(2);
    % v1 = R(1,1)^2;
    % v2 = R(1,2)^2;
    % v3 = R(1,3)^2;
    % v4 = R(1,4)^2;
    % v5 = R(1,5)^2;
    % v6 = R(1,6)^2;
    % v7 = R(1,7)^2;
    % v8 = R(1,8)^2;
    % u12 = R(1,1)*R(1,2);
    % u34 = R(1,3)*R(1,4);
    % u56 = R(1,5)*R(1,6);
    % u78 = R(1,7)*R(1,8);
    %xxx = [h b v1 v2 v3 v4 v5 v6 v7 v8 u12 u34 u56 u78]'; %ye
    
    %Ok, so this checks so that the signas extracted are "correct". If
    %we have complex solution, I just don't care what it does.
    if xxxsol(11) >0
        vv1 = sqrt(xxxsol(3)); % Here is something to fix. There are mirrored solutions
        vv2 = sqrt(xxxsol(4)); % One really have to check both +sqrt( ... ) and -sqrt( ...)
    else
        vv1 = sqrt(xxxsol(3));
        vv2 = -sqrt(xxxsol(4));
    end
    
    if xxxsol(12) >0
        vv3 = sqrt(xxxsol(5)); %
        vv4 = sqrt(xxxsol(6));
    else
        vv3 = sqrt(xxxsol(5));
        vv4 = -sqrt(xxxsol(6));
    end
    
    if xxxsol(13) > 0
        vv5 = sqrt(xxxsol(7));
        vv6 = sqrt(xxxsol(8));
    else
        vv5 = sqrt(xxxsol(7));
        vv6 = -sqrt(xxxsol(8));
    end
    
    if xxxsol(14) >0
        vv7 = sqrt(xxxsol(9));
        vv8 = sqrt(xxxsol(10));  %there are the values of x-coords for Receovers
    else
        vv7 = sqrt(xxxsol(9));
        vv8 = -sqrt(xxxsol(10));  %there are the values of x-coords for Receovers
    end
    
    
    
    %
    RRp = inv(ll)*Rtilde;
    SSp = (ll*(Stilde+bb))/(-2);
    % Also calculate S(3,1) and S(3,2)
    % using
    % d11*h - b^2/4 - h*S(3,1)^2 - h*R(1,1)'*R(1,1)
    % and
    % d1j(2)*h - (Stilde(2)*Stilde(2)/4 + 2*b*Stilde(2)/4 + h*S(3,2)^2 - h*S(3,1)^2)
    
    % S(3,1) = sqrt( d11 - b^2/(4*h - R(1,1)'*R(1,1) )
    S31 = sqrt(d11 - bb^2/(4*hh)  - vv1^2);
    % S(3,2) = sqrt( d1j(2) - (Stilde(2)*Stilde(2)/(4*h) + 2*b*Stilde(2)/(4*h) - S(3,1)^2)
    S32 = sqrt( d1j(2) - (Stilde(2)*Stilde(2)/(4*hh) + 2*bb*Stilde(2)/(4*hh) - S31^2) ); %OBS! HAVENT CHECKED THIS
    
    Rsol = [vv1 vv2 vv3 vv4 vv5 vv6 vv7 vv8;RRp;zeros(1,8)];
    Ssol = [0 0;SSp;S31 S32];
    
    %adding the different mirrored solutions
    %adding the different mirrored solutions
    
    %these for-loops are for gooing throug the 2^3 different mic rig
    %mirroring cases. Innermost we add the foru differnet sols for sender
    %mirrorings.
    %WARNING: THis adds quite a few solutions that are just completely
    %mirrored solutions in either y-z-plane or x-y plane. SO, getting 4
    %times as many solutions as there should be.
    
    %THE LAST MIC RIG  NEED NOT BE MIRRORED! THey will only get you solutions
    %that are complete mirrored in y-z-plane of the ones already added.
    for mm=0:1 %first rig
        if mm %mirror first rig
            Rsol(1,1) = -Rsol(1,1);
            Rsol(1,2) = -Rsol(1,2);
            
        end
        for nn=0:1 %second rig
            
            if nn %mirror second rig
                Rsol(1,3) = -Rsol(1,3);
                Rsol(1,4) = -Rsol(1,4);
            end
            
            for oo=0:1 %third rig
                if oo
                    Rsol(1,5) = -Rsol(1,5);
                    Rsol(1,6) = -Rsol(1,6);
                end
                
                %and add the solutions... 2 sols is added for each of the mirrored
                %senders
                sols.R{end+1}=Rsol;
                sols.S{end+1}=Ssol;
                
                Ssol(3,1)=-Ssol(3,1);
                sols.R{end+1}=Rsol;
                sols.S{end+1}=Ssol;
                %THESE TWO NEED NOT BE ADDED! As these are only mirrired
                %solutions of the two above, mirrored in the x-y-plane.
                %                 Ssol(3,2)=-Ssol(3,2);
                %                 sols.R{end+1}=Rsol;
                %                 sols.S{end+1}=Ssol;
                %
                %                 Ssol(3,1)=-Ssol(3,1);
                %                 sols.R{end+1}=Rsol;
                %                 sols.S{end+1}=Ssol;
                
            end
        end
    end
    
    
    
end
 sols.psols = psols;

end
%[err,errlist] = evaluate_solutions(sols,gt,settings);

%err
