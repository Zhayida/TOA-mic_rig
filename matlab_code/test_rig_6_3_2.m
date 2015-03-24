clear all

R = [zeros(2,1) rand(2,5);1+rand(1,6)];
R(2,2)=0;
S = [rand(2,3);zeros(1,3)];
%R = round(R*10);
%S = round(S*10);
Rp = R(1:2,:);
Sp = S(1:2,:);

d=toa_calc_d_from_xy(R,S);
d2 = d.^2;

d11 = d2(1,1);
d1j = d2(1,:)-d2(1,1);
di1 = d2(:,1)-d2(1,1);
dij = d2-repmat(d2(:,1),1,3)-repmat(d2(1,:),6,1)+d2(1,1);
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
% S(:,1)'*S(:,1) + R(3,1)'*R(3,1)
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
u12 = R(3,1)*R(3,2);
u34 = R(3,3)*R(3,4);
u56 = R(3,5)*R(3,6);


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

di1(2) - (Rtilde(:,2)'*H*Rtilde(:,2) -2*b'*Rtilde(:,2)/(-2) + v2 - v1)
di1(3) - (Rtilde(:,3)'*H*Rtilde(:,3) -2*b'*Rtilde(:,3)/(-2) + v3 - v1)
di1(4) - (Rtilde(:,4)'*H*Rtilde(:,4) -2*b'*Rtilde(:,4)/(-2) + v4 - v1)
di1(5) - (Rtilde(:,5)'*H*Rtilde(:,5) -2*b'*Rtilde(:,5)/(-2) + v5 - v1)
di1(6) - (Rtilde(:,6)'*H*Rtilde(:,6) -2*b'*Rtilde(:,6)/(-2) + v6 - v1)

clear A B
B(1,1) = di1(2);
B(2,1) = di1(3);
B(3,1) = di1(4);
B(4,1) = di1(5);
B(5,1) = di1(6);
A(1,1:14)=[Rtilde(1,2)^2 2*Rtilde(1,2)*Rtilde(2,2) Rtilde(2,2)^2 (-2*Rtilde(:,2)'/(-2)) -1 1 0 0 0 0 0 0 0];
A(2,1:14)=[Rtilde(1,3)^2 2*Rtilde(1,3)*Rtilde(2,3) Rtilde(2,3)^2 (-2*Rtilde(:,3)'/(-2)) -1 0 1 0 0 0 0 0 0];
A(3,1:14)=[Rtilde(1,4)^2 2*Rtilde(1,4)*Rtilde(2,4) Rtilde(2,4)^2 (-2*Rtilde(:,4)'/(-2)) -1 0 0 1 0 0 0 0 0];
A(4,1:14)=[Rtilde(1,5)^2 2*Rtilde(1,5)*Rtilde(2,5) Rtilde(2,5)^2 (-2*Rtilde(:,5)'/(-2)) -1 0 0 0 1 0 0 0 0];
A(5,1:14)=[Rtilde(1,6)^2 2*Rtilde(1,6)*Rtilde(2,6) Rtilde(2,6)^2 (-2*Rtilde(:,6)'/(-2)) -1 0 0 0 0 1 0 0 0];
xxx = [H(1,1) H(1,2) H(2,2) b(1) b(2) v1 v2 v3 v4 v5 v6 u12 u34 u56]';
A*xxx-B

%%%%%%%%%%%%%%%%
%% l^2
U12 = norm(R(:,2)-R(:,1));
U12^2
R(:,2)'*R(:,2)-2*R(:,2)'*R(:,1) + R(:,1)'*R(:,1)
Rtilde(:,2)'*H*Rtilde(:,2) + R(3,2)^2 + R(3,1)^2 - 2*R(3,1)*R(3,2)
B(6,1) = U12^2;
A(6,1:14)=[Rtilde(1,2)^2 2*Rtilde(1,2)*Rtilde(2,2) Rtilde(2,2)^2 0 0 1 1 0 0 0 0 -2 0 0];

%%%%%%%%%%%%%%%%
%% l^2
U34 = norm(R(:,4)-R(:,3));
U34^2
R(:,4)'*R(:,4)-2*R(:,3)'*R(:,4) + R(:,3)'*R(:,3)
tmp = Rtilde(:,4)*Rtilde(:,4)' -2*Rtilde(:,3)*Rtilde(:,4)' + Rtilde(:,3)*Rtilde(:,3)';
sum(sum(H.*tmp)) + R(3,4)^2 + R(3,3)^2 - 2*R(3,3)*R(3,4)
B(7,1) = U34^2;
A(7,1:14)=[tmp(1,1) (tmp(2,1)+tmp(1,2)) tmp(2,2) 0 0 0 0 1 1 0 0 0 -2 0];

%%%%%%%%%%%%%%%%
%% l^2
U56 = norm(R(:,6)-R(:,5));
U56^2
R(:,6)'*R(:,6)-2*R(:,5)'*R(:,6) + R(:,5)'*R(:,5)
tmp = Rtilde(:,6)*Rtilde(:,6)' -2*Rtilde(:,5)*Rtilde(:,6)' + Rtilde(:,5)*Rtilde(:,5)';
sum(sum(H.*tmp)) + R(3,6)^2 + R(3,5)^2 - 2*R(3,5)*R(3,6)
B(8,1) = U56^2;
A(8,1:14)=[tmp(1,1) (tmp(2,1)+tmp(1,2)) tmp(2,2) 0 0 0 0 0 0 1 1 0 0 -2];

[A*xxx B]

% OK. So now we have 8 linear equations and 6 non-linear equations.
% Three of the nonlinear are easy to express:
% v1*v2-u12^2 = 0
% v3*v4-u34^2 = 0
% v5*v6-u56^2 = 0
% But the three involving d11 and d1j involve inv(H) and are more
% difficult.
%
% v1*v2-u12^2 = 0
% v3*v4-u34^2 = 0
% v5*v6-u56^2 = 0
% b'*adj(H)*b/4+ det(H)*v1 - det(H)*d11
% Stilde(:,2)'*adj(H)*Stilde(:,2)/4 + 2*b'*adj(H)*Stilde(:,2)/4 - det(H)*d1j(2)
% Stilde(:,3)'*adj(H)*Stilde(:,3)/4 + 2*b'*adj(H)*Stilde(:,3)/4 - det(H)*d1j(3)



A1 = A(:,7:end);
A2 = A(:,1:6);
C1 = -A1\A2;
C2 = A1\B;
monvek(1,1)=multipol(1,[1;0;0;0;0;0]);
monvek(2,1)=multipol(1,[0;1;0;0;0;0]);
monvek(3,1)=multipol(1,[0;0;1;0;0;0]);
monvek(4,1)=multipol(1,[0;0;0;1;0;0]);
monvek(5,1)=multipol(1,[0;0;0;0;1;0]);
monvek(6,1)=multipol(1,[0;0;0;0;0;1]);
ett = multipol(1,[0;0;0;0;0;0]);
vu = C1*monvek+C2*ett;
HH = [monvek(1) monvek(2);monvek(2) monvek(3)];
bb = [monvek(4);monvek(5)];
vv1 = monvek(6);
eqs(1,1)=bb'*adjv(HH)*bb/4+ det(HH)*vv1 - detv(HH)*d11;
eqs(2,1)=Stilde(:,2)'*adjv(HH)*Stilde(:,2)/4 + 2*bb'*adjv(HH)*Stilde(:,2)/4 - det(HH)*d1j(2);
eqs(3,1)=Stilde(:,3)'*adjv(HH)*Stilde(:,3)/4 + 2*bb'*adjv(HH)*Stilde(:,3)/4 - det(HH)*d1j(3);
eqs(4,1)=monvek(6)*vu(1)-vu(6)^2;
eqs(5,1)=vu(2)*vu(3)-vu(7)^2;
eqs(6,1)=vu(4)*vu(5)-vu(8)^2;

% b'*adj(H)*b/4+ det(H)*v1 - det(H)*d11
% Stilde(:,2)'*adj(H)*Stilde(:,2)/4 + 2*b'*adj(H)*Stilde(:,2)/4 - det(H)*d1j(2)
% Stilde(:,3)'*adj(H)*Stilde(:,3)/4 + 2*b'*adj(H)*Stilde(:,3)/4 - det(H)*d1j(3)
% v1*v2-u12^2 = 0
% v3*v4-u34^2 = 0
% v5*v6-u56^2 = 0

gt = xxx(1:6);
[evaluate(vu,gt) xxx(7:end)]
evaluate(eqs,gt)

eqs2m2(eqs*1000000,'testrig63.m2');


if 1
    %%
    
    [mm,nn]  = get_degree(eqs);
    
    NN = 4;
    MM = 6;
    mmax = [2 3 3 NN NN 1]';
    tmax = [MM MM MM MM MM MM];

    
    C = polynomials2matrix(eqs);
    eqs2 = [];
    for i = 1:length(eqs(1:end));
        multimons{i} = monvec(create_upto_degree(mmax-mm{i},tmax(i)-nn(i)));
        eqs(i) = eqs(i)/norm(C(i,:));
        eqs2 = [eqs2(:);eqs(i)*multimons{i}(:)];
    end
    % eqs2 = generate_equations(eqs,4);
    
    settings.dim = 16;
    settings.basis_size = 30;
    % settings.permissible = B2;
    settings.action_variable = 1;
    [sols stats] = polysolve(eqs2, settings);
    err = evaluate_solutions(sols,gt,settings)
    stats
end
