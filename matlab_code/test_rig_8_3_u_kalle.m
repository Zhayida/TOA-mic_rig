clearvars

integerproblem = 1;

%addpath('C:\Users\simonb\Documents\antennaResection\matlab\bundleadjusters')
%% Generate Sythetic data

%generating integer ground truth, for integer measurements.
mMics = 8;
nSounds = 3;
rigDist = 150;
factor=20;

R = factor*[zeros(2,1) rand(2,7); 1 + rand(1,mMics)];
R(2,2) = 0;
R = round(R);

if 0,
    for kk = 1 : (size(R,2)/2)
        phi = 2*pi*rand;
        Rot = [cos(phi) sin(phi); -sin(phi) cos(phi)];
        R(:,kk*2) = R(:,kk*2-1) + [Rot*rigDist*[1; 0] ; 0] ;
    end
else
    for kk = 1 : (size(R,2)/2)
        deltaR = [3;4;1];
        R(:,kk*2) = R(:,kk*2-1) + deltaR ;
    end   
end;

%%Sanity check on receiver ditances = rigDist
norm(R(:,1)-R(:,2))
norm(R(:,3)-R(:,4))
norm(R(:,5)-R(:,6))
norm(R(:,7)-R(:,8))


S = factor*[rand(2,3); zeros(1,3)];
S = round(S);
% R = round(R*10);
% S = round(S*10);

if 0,
    
end

%% Prameterization
Rp = R(1:2,:); %Receiver projection ground truth
Sp = S(1:2,:); %Sender projection grounf truth

d = toa_calc_d_from_xy(R,S);

%dRound=round(d);
%[I,J,D]=find(dRound);
%[xopt,yopt]=bundleToaMicRigsUncalib(D,I,J,R,S);
%dTest = toa_calc_d_from_xy(xopt,yopt);

d2 = d.^2;

d11 = d2(1,1);
d1j = d2(1,:)-d2(1,1);
di1 = d2(:,1)-d2(1,1);
dij = d2-repmat(d2(:,1),1,3)-repmat(d2(1,:),8,1)+d2(1,1);
Rtilde = [zeros(2,1) dij(2:end,2:3)'];
Stilde = [zeros(2,1) eye(2)];

L = (Rtilde/Rp)';
b1 = mean( (L*Stilde/(-2)-Sp)' )';
b = 2*inv(L)*b1;
Rp
inv(L')*Rtilde
[Rp - inv(L')*Rtilde]
Sp
(L*(Stilde+repmat(b,1,3)))/(-2)

H = inv(L'*L)



%%%%%%%%%%%%%%%%
%%% d11
d11
S(:,1)'*S(:,1) + R(3,1)'*R(3,1)
% b'*inv(H)*b/4+ R(3,1)'*R(3,1)
% b'*adjv(H)*b/4+ det(H)*v1 - det(H)*d11

%%%%%%%%%%%%%%%%
%%% d1j
d1j(2)
S(:,2)'*S(:,2)-S(:,1)'*S(:,1)
Stilde(:,2)'*inv(H)*Stilde(:,2)/4 + 2*b'*inv(H)*Stilde(:,2)/4
Stilde(:,2)'*adjv(H)*Stilde(:,2)/4 + 2*b'*adjv(H)*Stilde(:,2)/4 - det(H)*d1j(2)
%d1j(2)*h - (Stilde(2)*Stilde(2)/4 + 2*b*Stilde(2)/4 + h*S(3,2)^2 - h*S(3,1)^2)
d1j(3)
S(:,3)'*S(:,3)-S(:,1)'*S(:,1)
Stilde(:,3)'*inv(H)*Stilde(:,3)/4 + 2*b'*inv(H)*Stilde(:,3)/4
Stilde(:,3)'*adjv(H)*Stilde(:,3)/4 + 2*b'*adjv(H)*Stilde(:,3)/4 - det(H)*d1j(3)
%d1j(2)*h - (Stilde(2)*Stilde(2)/4 + 2*b*Stilde(2)/4 + h*S(3,2)^2 - h*S(3,1)^2)

v1 = R(3,1)^2;
v2 = R(3,2)^2;
v3 = R(3,3)^2;
v4 = R(3,4)^2;
v5 = R(3,5)^2;
v6 = R(3,6)^2;
v7 = R(3,7)^2;
v8 = R(3,8)^2;
u12 = R(3,1)*R(3,2);
u34 = R(3,3)*R(3,4);
u56 = R(3,5)*R(3,6);
u78 = R(3,7)*R(3,8);

%%%%%%%%%%%%%%%%
%%% di1
di1(2)
R(:,2)'*R(:,2) - 2*R(:,2)'*S(:,1) - R(:,1)'*R(:,1)
Rtilde(:,2)'*H*Rtilde(:,2) -2*b'*Rtilde(:,2)/(-2) + R(3,2)^2 - R(3,1)^2
di1(2) - (Rtilde(:,2)'*H*Rtilde(:,2) -2*b'*Rtilde(:,2)/(-2) + R(3,2)^2 - R(3,1)^2)
di1(3) - (Rtilde(:,3)'*H*Rtilde(:,3) -2*b'*Rtilde(:,3)/(-2) + R(3,3)^2 - R(3,1)^2)
di1(4) - (Rtilde(:,4)'*H*Rtilde(:,4) -2*b'*Rtilde(:,4)/(-2) + R(3,4)^2 - R(3,1)^2)
di1(5) - (Rtilde(:,5)'*H*Rtilde(:,5) -2*b'*Rtilde(:,5)/(-2) + R(3,5)^2 - R(3,1)^2)
di1(6) - (Rtilde(:,6)'*H*Rtilde(:,6) -2*b'*Rtilde(:,6)/(-2) + R(3,6)^2 - R(3,1)^2)
di1(7) - (Rtilde(:,7)'*H*Rtilde(:,7) -2*b'*Rtilde(:,7)/(-2) + R(3,7)^2 - R(3,1)^2)
di1(8) - (Rtilde(:,8)'*H*Rtilde(:,8) -2*b'*Rtilde(:,8)/(-2) + R(3,8)^2 - R(3,1)^2)

di1(2) - (Rtilde(:,2)'*H*Rtilde(:,2) -2*b'*Rtilde(:,2)/(-2) + v2 - v1)
di1(3) - (Rtilde(:,3)'*H*Rtilde(:,3) -2*b'*Rtilde(:,3)/(-2) + v3 - v1)
di1(4) - (Rtilde(:,4)'*H*Rtilde(:,4) -2*b'*Rtilde(:,4)/(-2) + v4 - v1)
di1(5) - (Rtilde(:,5)'*H*Rtilde(:,5) -2*b'*Rtilde(:,5)/(-2) + v5 - v1)
di1(6) - (Rtilde(:,6)'*H*Rtilde(:,6) -2*b'*Rtilde(:,6)/(-2) + v6 - v1)
di1(7) - (Rtilde(:,7)'*H*Rtilde(:,7) -2*b'*Rtilde(:,7)/(-2) + v7 - v1)
di1(8) - (Rtilde(:,8)'*H*Rtilde(:,8) -2*b'*Rtilde(:,8)/(-2) + v8 - v1)

clear A B
B(1,1) = di1(2);
B(2,1) = di1(3);
B(3,1) = di1(4);
B(4,1) = di1(5);
B(5,1) = di1(6);
B(6,1) = di1(7);
B(7,1) = di1(8);

A(1,1:17)=[Rtilde(1,2)^2 2*Rtilde(1,2)*Rtilde(2,2) Rtilde(2,2)^2 (-2*Rtilde(:,2)'/(-2)) -1 1 0 0 0 0 0 0 0 0 0 0];
A(2,1:17)=[Rtilde(1,3)^2 2*Rtilde(1,3)*Rtilde(2,3) Rtilde(2,3)^2 (-2*Rtilde(:,3)'/(-2)) -1 0 1 0 0 0 0 0 0 0 0 0];
A(3,1:17)=[Rtilde(1,4)^2 2*Rtilde(1,4)*Rtilde(2,4) Rtilde(2,4)^2 (-2*Rtilde(:,4)'/(-2)) -1 0 0 1 0 0 0 0 0 0 0 0];
A(4,1:17)=[Rtilde(1,5)^2 2*Rtilde(1,5)*Rtilde(2,5) Rtilde(2,5)^2 (-2*Rtilde(:,5)'/(-2)) -1 0 0 0 1 0 0 0 0 0 0 0];
A(5,1:17)=[Rtilde(1,6)^2 2*Rtilde(1,6)*Rtilde(2,6) Rtilde(2,6)^2 (-2*Rtilde(:,6)'/(-2)) -1 0 0 0 0 1 0 0 0 0 0 0];
A(6,1:17)=[Rtilde(1,7)^2 2*Rtilde(1,7)*Rtilde(2,7) Rtilde(2,7)^2 (-2*Rtilde(:,7)'/(-2)) -1 0 0 0 0 0 1 0 0 0 0 0];
A(7,1:17)=[Rtilde(1,8)^2 2*Rtilde(1,8)*Rtilde(2,8) Rtilde(2,8)^2 (-2*Rtilde(:,8)'/(-2)) -1 0 0 0 0 0 0 1 0 0 0 0];

xxx = [H(1,1) H(1,2) H(2,2) b(1) b(2) v1 v2 v3 v4 v5 v6 v7 v8 u12 u34 u56 u78]'
%% l^2
% we have ||r1-r2||=||r3-r4||=||r5-r6||=||r7-r8||= l^2

% use equation ||r1-r2||-||r3-r4||=0
norm(R(:,2)-R(:,1))- norm(R(:,4)-R(:,3))

%groud truth test
u12 = R(:,2)'*R(:,2) - 2*R(:,2)'*R(:,1) + R(:,1)'*R(:,1);
u34 = R(:,3)'*R(:,3) - 2*R(:,3)'*R(:,4) + R(:,4)'*R(:,4);
u12 - u34

% Parameterization of ground truth
uu12 = Rtilde(:,2)'*H*Rtilde(:,2) + R(3,2)^2 + R(3,1)^2 - 2*R(3,1)*R(3,2);
tmp = Rtilde(:,4)*Rtilde(:,4)' -2*Rtilde(:,3)*Rtilde(:,4)' + Rtilde(:,3)*Rtilde(:,3)';
uu34 = sum(sum(H.*tmp)) + R(3,4)^2 + R(3,3)^2 - 2*R(3,3)*R(3,4);
uu12 - uu34

B(8,1) = 0; %we take the two rig distances minus each other
a1 = Rtilde(1,2)^2 - tmp(1,1);
a2 = 2*Rtilde(1,2)*Rtilde(2,2) - (tmp(2,1)+tmp(1,2));
a3 = Rtilde(2,2)^2 - tmp(2,2);
A(8,1:17)=[a1 a2 a3 0 0 1 1 -1 -1 0 0 0 0 -2 2 0 0];

%% l^2
% use equation ||r1-r2||-||r5-r6||=0
norm(R(:,2)-R(:,1))- norm(R(:,6)-R(:,5))

%groud truth test
u56 = R(:,5)'*R(:,5) - 2*R(:,5)'*R(:,6) + R(:,6)'*R(:,6);
u12 - u56

% Parameterization of ground truth
tmp = Rtilde(:,6)*Rtilde(:,6)' -2*Rtilde(:,5)*Rtilde(:,6)' + Rtilde(:,5)*Rtilde(:,5)';
uu56 = sum(sum(H.*tmp)) + R(3,6)^2 + R(3,5)^2 - 2*R(3,5)*R(3,6)
uu12 - uu56

B(9,1) = 0; %we take the two rig distances minus each other
a1 = Rtilde(1,2)^2 - tmp(1,1);
a2 = 2*Rtilde(1,2)*Rtilde(2,2) - (tmp(2,1)+tmp(1,2));
a3 = Rtilde(2,2)^2 - tmp(2,2);
A(9,1:17)=[a1 a2 a3 0 0 1 1 0 0 -1 -1 0 0 -2 0 2 0];
%% l^2
% use equation ||r1-r2||-||r7-r8||=0
norm(R(:,2)-R(:,1))- norm(R(:,8)-R(:,7))

%groud truth test
u56 = R(:,7)'*R(:,7) - 2*R(:,7)'*R(:,8) + R(:,8)'*R(:,8);
u12 - u56

% Parameterization of ground truth
tmp = Rtilde(:,8)*Rtilde(:,8)' -2*Rtilde(:,7)*Rtilde(:,8)' + Rtilde(:,7)*Rtilde(:,7)';
uu78 = sum(sum(H.*tmp)) + R(3,8)^2 + R(3,7)^2 - 2*R(3,7)*R(3,8)
uu12 - uu78

B(10,1) = 0; %we take the two rig distances minus each other
a1 = Rtilde(1,2)^2 - tmp(1,1);
a2 = 2*Rtilde(1,2)*Rtilde(2,2) - (tmp(2,1)+tmp(1,2));
a3 = Rtilde(2,2)^2 - tmp(2,2);
A(10,1:17)=[a1 a2 a3 0 0 1 1 0 0 0 0 -1 -1 -2 0 0 2];
%%
norm(A*xxx-B)
A = round(A);
b = round(b);

%ids = [7:16];
ids = [7:13 15:17];
rids = setdiff(1:17,ids);
A1 = A(:,ids);
A2 = A(:,rids);
C1 = -A1\A2;
C2 = A1\B;
E = eye(7);

for i = 1:7
    monvek(i,1)=multipol(1,E(:,i));
end

ett = multipol(1,[0;0;0;0;0;0;0]);
vu = C1*monvek+C2*ett;
HH = [monvek(1) monvek(2);monvek(2) monvek(3)];
bb = [monvek(4);monvek(5)];
vv1 = monvek(6);
u78 = monvek(7);
eqs(1,1)=bb'*adjv(HH)*bb/4+ det(HH)*vv1 - detv(HH)*d11;
eqs(2,1)=Stilde(:,2)'*adjv(HH)*Stilde(:,2)/4 + 2*bb'*adjv(HH)*Stilde(:,2)/4 - det(HH)*d1j(2);
eqs(3,1)=Stilde(:,3)'*adjv(HH)*Stilde(:,3)/4 + 2*bb'*adjv(HH)*Stilde(:,3)/4 - det(HH)*d1j(3);
eqs(4,1)=vv1*vu(1)-u78^2;
eqs(5,1)=vu(2)*vu(3)-vu(8)^2;
eqs(6,1)=vu(4)*vu(5)-vu(9)^2;
eqs(7,1)=vu(6)*vu(7)-vu(10)^2;


% b'*adj(H)*b/4+ det(H)*v1 - det(H)*d11
% Stilde(:,2)'*adj(H)*Stilde(:,2)/4 + 2*b'*adj(H)*Stilde(:,2)/4 - det(H)*d1j(2)
% Stilde(:,3)'*adj(H)*Stilde(:,3)/4 + 2*b'*adj(H)*Stilde(:,3)/4 - det(H)*d1j(3)
% v1*v2-u12^2 = 0
% v3*v4-u34^2 = 0
% v5*v6-u56^2 = 0

gt = xxx(rids);
% [evaluate(vu,gt) xxx(7:end)]
evaluate(eqs,gt)

eqs2m2(eqs*4,'testrig83kalle.m2');
%testNbrSols_M2(eqs);
eval(['!m2<','testrig83kalle.m2'])
% 29 solutions
% Basis of quotient ring according to macauley
% 1 x1 x1x2 x1x3 x1x5 x1x7 x2 x2x3 x2x5 x2x6 x2x7 x3 x3^2 x3x5 x3x6 x3x7 x4 x4x5 x4x6 x4x7 x5 x5^2 x5x6 x5x7 x6 x6^2 x6x7 x7 x7^2 
% leading term of groebner basis
lta = ['x1x6 x4^2 x3x4 x2x4 x1x4 x2^2 x1^2 x7^3 x6x7^2 x5x7^2 x4x7^2 x3x7^2 x2x7^2 x1x7^2 x6^2x7 x5x6x7 ' ...
'x4x6x7 x3x6x7 x2x6x7 x5^2x7 x4x5x7 x3x5x7 x2x5x7 x1x5x7 x3^2x7 x2x3x7 x1x3x7 x1x2x7 x6^3 x5x6^2 x4x6^2 ' ...
'x3x6^2 x2x6^2 x5^2x6 x4x5x6 x3x5x6 x2x5x6 x3^2x6 x2x3x6 x5^3 x4x5^2 x3x5^2 x2x5^2 x1x5^2 x3^2x5 x2x3x5 ' ...
'x1x3x5 x1x2x5 x3^3 x2x3^2 x1x3^2 x1x2x3 '];

% return
%%
if 1
    [mm,nn]  = get_degree(eqs);
    
    NN = 3;
    MM = 7;
    mmax = [NN NN NN NN NN NN 3]';
    tmax = [MM MM MM MM MM MM MM];
    
    
    C = polynomials2matrix(eqs);
    eqs2 = [];
    for i = 1:length(eqs(1:end));
        multimons{i} = monvec(create_upto_degree(mmax-mm{i},tmax(i)-nn(i)));
        eqs(i) = eqs(i)/norm(C(i,:));
        eqs2 = [eqs2(:);eqs(i)*multimons{i}(:)];
    end
    % eqs2 = generate_equations(eqs,4);
    
    settings.dim = 80;
    settings.basis_size = 130;
    % settings.permissible = B2;
    settings.action_variable = 1;
    [sols stats] = polysolve(eqs2, settings);
    err = evaluate_solutions(sols,gt,settings)
    stats
end