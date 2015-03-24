clearvars
close all
if ispc
addpath .
addpath .

end
%%
R = [zeros(2,1) rand(2,5);zeros(1,6)];
R(1,1)=rand(1,1);
S = [0 0;rand(1,1)*0.2 1+rand(1,1);rand(1,2)];
Rp = R(2,:);
Sp = S(2,:);

d=toa_calc_d_from_xy(R,S);
d2 = d.^2;

d11 = d2(1,1);
d1j = d2(1,:)-d2(1,1);
di1 = d2(:,1)-d2(1,1);
dij = d2-repmat(d2(:,1),1,2)-repmat(d2(1,:),6,1)+d2(1,1);
Rtilde = [0 dij(2:end,2)'];
Stilde = [0 1];

l = Rtilde/Rp;
b = mean(-2*Sp/(l))-mean(Stilde) %So s_(2,1)=l*b/-2 and s_(2,2)=l*(Stilde(2)+b)/-2
Rp
inv(l)*Rtilde
Sp 
(l*(Stilde+b))/(-2)
h = 1/(l^2)


%%%%%%%%%%%%%%%%
%%% d11
d11
S(:,1)'*S(:,1)+ R(1,1)'*R(1,1)
b*l*l*b/4+S(3,1)^2 + R(1,1)'*R(1,1)
d11*h - b^2/4 - h*S(3,1)^2 - h*R(1,1)'*R(1,1)

%%%%%%%%%%%%%%%%
%%% d1j
d1j(2)
S(:,2)'*S(:,2)-S(:,1)'*S(:,1)
Stilde(2)*l*l*Stilde(2)/4 + 2*b*l*l*Stilde(2)/4 + S(3,2)^2 - S(3,1)^2
d1j(2)*h - (Stilde(2)*Stilde(2)/4 + 2*b*Stilde(2)/4 + h*S(3,2)^2 - h*S(3,1)^2)

v1 = R(1,1)^2;
v2 = R(1,2)^2;
v3 = R(1,3)^2;
v4 = R(1,4)^2;
v5 = R(1,5)^2;
v6 = R(1,6)^2;
u12 = R(1,1)*R(1,2);
u34 = R(1,3)*R(1,4);
u56 = R(1,5)*R(1,6);


%%%%%%%%%%%%%%%%
%%% di1
di1(2)
R(:,2)'*R(:,2) - 2*R(:,2)'*S(:,1) - R(:,1)'*R(:,1)
Rtilde(2)*inv(l*l)*Rtilde(2) -2*b'*Rtilde(2)/(-2) + R(1,2)^2 - R(1,1)^2
di1(2) - (Rtilde(2)*inv(l*l)*Rtilde(2) -2*b'*Rtilde(2)/(-2) + R(1,2)^2 - R(1,1)^2)
di1(3) - (Rtilde(3)*inv(l*l)*Rtilde(3) -2*b'*Rtilde(3)/(-2) + R(1,3)^2 - R(1,1)^2)
di1(4) - (Rtilde(4)*inv(l*l)*Rtilde(4) -2*b'*Rtilde(4)/(-2) + R(1,4)^2 - R(1,1)^2)
di1(5) - (Rtilde(5)*inv(l*l)*Rtilde(5) -2*b'*Rtilde(5)/(-2) + R(1,5)^2 - R(1,1)^2)
di1(6) - (Rtilde(6)*inv(l*l)*Rtilde(6) -2*b'*Rtilde(6)/(-2) + R(1,6)^2 - R(1,1)^2)

di1(2) - (Rtilde(2)*h*Rtilde(2) -2*b'*Rtilde(2)/(-2) + v2 - v1)
di1(3) - (Rtilde(3)*h*Rtilde(3) -2*b'*Rtilde(3)/(-2) + v3 - v1)
di1(4) - (Rtilde(4)*h*Rtilde(4) -2*b'*Rtilde(4)/(-2) + v4 - v1)
di1(5) - (Rtilde(5)*h*Rtilde(5) -2*b'*Rtilde(5)/(-2) + v5 - v1)
di1(6) - (Rtilde(6)*h*Rtilde(6) -2*b'*Rtilde(6)/(-2) + v6 - v1)

clear A B
B(1,1) = di1(2);
B(2,1) = di1(3);
B(3,1) = di1(4);
B(4,1) = di1(5);
B(5,1) = di1(6); %Unknowns are in order of ( h  b v1 v2 v3 v4 v5 v6
A(1,1:11)=[Rtilde(2)^2 (-2*Rtilde(2)/(-2)) -1 1 0 0 0 0 0 0 0];
A(2,1:11)=[Rtilde(3)^2 (-2*Rtilde(3)/(-2)) -1 0 1 0 0 0 0 0 0];
A(3,1:11)=[Rtilde(4)^2 (-2*Rtilde(4)/(-2)) -1 0 0 1 0 0 0 0 0];
A(4,1:11)=[Rtilde(5)^2 (-2*Rtilde(5)/(-2)) -1 0 0 0 1 0 0 0 0];
A(5,1:11)=[Rtilde(6)^2 (-2*Rtilde(6)/(-2)) -1 0 0 0 0 1 0 0 0];
xxx = [h b v1 v2 v3 v4 v5 v6 u12 u34 u56]';


%%%%%%%%%%%%%%%%
%% ridDist^2 -equations
U12 = norm(R(:,2)-R(:,1));
U12^2
R(:,2)'*R(:,2)-2*R(:,2)'*R(:,1) + R(:,1)'*R(:,1)
Rtilde(2)*inv(l*l)*Rtilde(2) + R(1,2)^2 + R(1,1)^2 - 2*R(1,1)*R(1,2)
B(6,1) = U12^2;
A(6,1:11)=[Rtilde(2)^2 0 1 1 0 0 0 0 -2 0 0];

%%%%%%%%%%%%%%%%
%% l^2
U34 = norm(R(:,4)-R(:,3));
U34^2
R(:,4)'*R(:,4)-2*R(:,4)'*R(:,3) + R(:,3)'*R(:,3)
Rtilde(4)*inv(l*l)*Rtilde(4) - 2*Rtilde(3)*inv(l*l)*Rtilde(4) + Rtilde(3)*inv(l*l)*Rtilde(3) + R(1,4)^2 + R(1,3)^2 - 2*R(1,4)*R(1,3)
B(7,1) = U34^2;
A(7,1:11)=[(Rtilde(4)^2-2*Rtilde(3)*Rtilde(4) + Rtilde(3)^2) 0 0 0 1 1 0 0 0 -2 0];

%%%%%%%%%%%%%%%%
%% l^2
U56 = norm(R(:,6)-R(:,5));
U56^2
R(:,6)'*R(:,6)-2*R(:,6)'*R(:,5) + R(:,5)'*R(:,5)
Rtilde(6)*inv(l*l)*Rtilde(6) - 2*Rtilde(6)*inv(l*l)*Rtilde(5) + Rtilde(5)*inv(l*l)*Rtilde(5) + R(1,6)^2 + R(1,5)^2 - 2*R(1,6)*R(1,5)
B(8,1) = U56^2;
A(8,1:11)=[(Rtilde(5)^2-2*Rtilde(5)*Rtilde(6) + Rtilde(6)^2) 0 0 0 0 0 1 1 0 0 -2];

[A*xxx B]

A1 = A(:,4:end);
A2 = A(:,1:3);
C1 = -A1\A2;
C2 = A1\B;
monvek(1,1)=multipol(1,[1;0;0]);
monvek(2,1)=multipol(1,[0;1;0]);
monvek(3,1)=multipol(1,[0;0;1]);
ett = multipol(1,[0;0;0]);
vu = C1*monvek+C2*ett;
eqs(1,1)=monvek(3)*vu(1)-vu(6)^2;
eqs(2,1)=vu(2)*vu(3)-vu(7)^2;
eqs(3,1)=vu(4)*vu(5)-vu(8)^2;
gt = xxx(1:3);
norm([evaluate(vu,gt) - xxx(4:end)])/norm(xxx(4:end))
evaluate(eqs,gt)

eqs2 = expandtodegree(eqs,[2 2 2],4);
settings.dim = 2;
settings.basis_size = 8;
[sols stats] = polysolve(eqs2, settings)

% So now we have basically got two possible solutions
% Maybe only one is feasible.

% for kk = 1:size(sols,2);
kk = 2; % Hardcoded. 
onesol = sols(:,kk);
xxxsol = [onesol;evaluate(vu,onesol)];
% xxx = [h b v1 v2 v3 v4 v5 v6 u12 u34 u56]';
hh = xxxsol(1);
ll = 1/sqrt(hh);
bb = xxxsol(2);
vv1 = sqrt(xxxsol(3)); % Here is something to fix. There are mirrored solutions
vv2 = sqrt(xxxsol(4)); % One really have to check both +sqrt( ... ) and -sqrt( ...)
vv3 = sqrt(xxxsol(5)); % and choose appropriate solutions.
vv4 = sqrt(xxxsol(6));
vv5 = sqrt(xxxsol(7));
vv6 = sqrt(xxxsol(8));
%
RRp = inv(ll)*Rtilde;
SSp = (ll*(Stilde+bb))/(-2);
% Also calculate S(3,1) and S(3,2)
% using
% d11*h - b^2/4 - h*S(3,1)^2 - h*R(1,1)'*R(1,1)
% and
% d1j(2)*h - (Stilde(2)*Stilde(2)/4 + 2*b*Stilde(2)/4 + h*S(3,2)^2 - h*S(3,1)^2)

% S(3,1) = sqrt( d11 - b^2/(4*h - R(1,1)'*R(1,1) )
S31 = sqrt(d11 - bb^2/(4*hh)  - vv1^2)
% S(3,2) = sqrt( d1j(2) - (Stilde(2)*Stilde(2)/(4*h) + 2*b*Stilde(2)/(4*h) - S(3,1)^2)
S32 = sqrt( d1j(2) - (Stilde(2)*Stilde(2)/(4*hh) + 2*bb*Stilde(2)/(4*hh) - S31^2) )


Rsol = [vv1 vv2 vv3 vv4 vv5 vv6;RRp;zeros(1,6)];
Ssol = [0 0;SSp;S31 S32];
addpath('C:\Users\simonb\Documents\antennaResection\matlab\sharedRoutines')
us = [U12 U34 U56];
[sols]=toa_rig_6_2(d,us);

