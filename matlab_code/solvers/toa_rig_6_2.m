function [sols]=toa_rig_6_2(d,rigDist);
%WARNING! THis code assumes that the solutions are the positive mirrorings.

% Square all distances
d2 = d.^2;

% Form the four types of distance squared differences
d11 = d2(1,1);                                          % Type A
d1j = d2(1,:)-d2(1,1);                                  % Type B
di1 = d2(:,1)-d2(1,1);                                  % Type C
dij = d2-repmat(d2(:,1),1,2)-repmat(d2(1,:),6,1)+d2(1,1); % Type D

% Use dij to establish Rtilde and Stilde (up to unknown l and b)
% i e Type D equations
Rtilde = [0 dij(2:end,2)'];
Stilde = [0 1];

% Get linear equations in x = [h b v1 v2 v3 v4 v5 v6 u12 u34 u56]';
%clear A B

% v1 = R(1,1)^2;
% v2 = R(1,2)^2;
% v3 = R(1,3)^2;
% v4 = R(1,4)^2;
% v5 = R(1,5)^2;
% v6 = R(1,6)^2;
% u12 = R(1,1)*R(1,2);
% u34 = R(1,3)*R(1,4);
% u56 = R(1,5)*R(1,6);

% Equations from di1 (Type C)
B(1,1) = di1(2);
B(2,1) = di1(3);
B(3,1) = di1(4);
B(4,1) = di1(5);
B(5,1) = di1(6);
A(1,1:11)=[Rtilde(2)^2 (-2*Rtilde(2)/(-2)) -1 1 0 0 0 0 0 0 0];
A(2,1:11)=[Rtilde(3)^2 (-2*Rtilde(3)/(-2)) -1 0 1 0 0 0 0 0 0];
A(3,1:11)=[Rtilde(4)^2 (-2*Rtilde(4)/(-2)) -1 0 0 1 0 0 0 0 0];
A(4,1:11)=[Rtilde(5)^2 (-2*Rtilde(5)/(-2)) -1 0 0 0 1 0 0 0 0];
A(5,1:11)=[Rtilde(6)^2 (-2*Rtilde(6)/(-2)) -1 0 0 0 0 1 0 0 0];

% Equations from rig constraints
U12 = rigDist(1);
U34 = rigDist(2);
U56 = rigDist(3);
B(6,1) = U12^2;
A(6,1:11)=[Rtilde(2)^2 0 1 1 0 0 0 0 -2 0 0];
B(7,1) = U34^2;
A(7,1:11)=[(Rtilde(4)^2-2*Rtilde(3)*Rtilde(4) + Rtilde(3)^2) 0 0 0 1 1 0 0 0 -2 0];
B(8,1) = U56^2;
A(8,1:11)=[(Rtilde(5)^2-2*Rtilde(5)*Rtilde(6) + Rtilde(6)^2) 0 0 0 0 0 1 1 0 0 -2];

% Express [v2 v3 v4 v5 v6 u12 u34 u56] linearly in terms of [h b v1]
A1 = A(:,4:end);
A2 = A(:,1:3);
C1 = -A1\A2;
C2 = A1\B;
monvek(1,1)=multipol(1,[1;0;0]);
monvek(2,1)=multipol(1,[0;1;0]);
monvek(3,1)=multipol(1,[0;0;1]);
ett = multipol(1,[0;0;0]);
vu = C1*monvek+C2*ett;

% Form quadratic equations u12^2 = v1*v2, u34^=v3*v4, u56^2 = v5*v6
eqs(1,1)=monvek(3)*vu(1)-vu(6)^2;
eqs(2,1)=vu(2)*vu(3)-vu(7)^2;
eqs(3,1)=vu(4)*vu(5)-vu(8)^2;

% Solve equations
eqs2 = expandtodegree(eqs,[2 2 2],4);
settings.dim = 2;
settings.basis_size = 8;
settings.action_variable = 1;
%keyboard;
[psols stats] = polysolve(eqs2, settings);

% So now we have basically got two possible solutions
% Maybe only one is feasible.
    sols=struct;
    sols.R=cell(0);
    sols.S=cell(0);
for kk = 1:size(psols,2);
    % For each solution try to backtrack and calculate all
    % unknown parameters
    onesol = psols(:,kk);
    xxxsol = [onesol;evaluate(vu,onesol)];
    % xxx = [h b v1 v2 v3 v4 v5 v6 u12 u34 u56]';
    hh = xxxsol(1);
    ll = 1/sqrt(hh);
    bb = xxxsol(2);
    %     vv1 = sqrt(xxxsol(3)); % Here is something to fix. There are mirrored solutions
    %     vv2 = sqrt(xxxsol(4)); % One really have to check both +sqrt( ... ) and -sqrt( ...)
    %     vv3 = sqrt(xxxsol(5)); % and choose appropriate solutions.
    %     vv4 = sqrt(xxxsol(6));
    %     vv5 = sqrt(xxxsol(7));
    %     vv6 = sqrt(xxxsol(8));
    
    
    %these if blocks test the signs of the paris of R(1,i) for a rig such
    %that
    % v1 = R(1,1)^2;
    % v2 = R(1,2)^2;
    % v3 = R(1,3)^2;
    % v4 = R(1,4)^2;
    % v5 = R(1,5)^2;
    % v6 = R(1,6)^2;
    % u12 = R(1,1)*R(1,2);
    % u34 = R(1,3)*R(1,4);
    % u56 = R(1,5)*R(1,6);

    if xxxsol(9) >0
        vv1 = sqrt(xxxsol(3)); % Here is something to fix. There are mirrored solutions
        vv2 = sqrt(xxxsol(4)); % One really have to check both +sqrt( ... ) and -sqrt( ...). FIXED BELOW
    else
        vv1 = sqrt(xxxsol(3));
        vv2 = -sqrt(xxxsol(4));
    end
    
    if xxxsol(10) >0
        vv3 = sqrt(xxxsol(5)); %
        vv4 = sqrt(xxxsol(6));
    else
        vv3 = sqrt(xxxsol(5));
        vv4 = -sqrt(xxxsol(6));
    end
    
    if xxxsol(11) > 0
        vv5 = sqrt(xxxsol(7));
        vv6 = sqrt(xxxsol(8));
    else
        vv5 = sqrt(xxxsol(7));
        vv6 = -sqrt(xxxsol(8));
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
    S32 = sqrt( d1j(2) - (Stilde(2)*Stilde(2)/(4*hh) + 2*bb*Stilde(2)/(4*hh) - S31^2) );
    
    Rsol = [vv1 vv2 vv3 vv4 vv5 vv6;RRp;zeros(1,6)];
    Ssol = [0 0;SSp;S31 S32];
    
    
    %adding the different mirrored solutions
    
    %these for-loops are for gooing throug the 2^3 different mic rig
    %mirroring cases. Innermost we add the foru differnet sols for sender
    %mirrorings.
    %WARNING: THis adds quite a few solutions that are just completely
    %mirrored solutions in either y-z-plane or x-y plane. SO, getting 4
    %times as many solutions as there should be.

             %THIS FIRST MIC RIG  NEED NOT BE MIRRIRED! THey will only get you solutions
            %that are completely mirrored in y-z-plane
    %     for mm=0:1 %first rig
%         if mm %mirror first rig
%             Rsol(1,1) = -Rsol(1,1);
%             Rsol(1,2) = -Rsol(1,2);
%             
       %    end
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
                
                %and add the solutions... 4 sols is added for each of the mirrored
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
    %end
    
    
    
     
%     sols.R{kk}=Rsol;
%     sols.S{kk}=Ssol;
    
end
sols.psols = psols;
end

