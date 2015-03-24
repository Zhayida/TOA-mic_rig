clearvars
close all
if ispc
addpath C:\Users\simonb\Documents\multipol
addpath .

end

%% Constants
findIntegerMeasurements=1;

%%
nMics=8;
nSounds=2;
%Creating groud truth
rigDist=150;
factor=1000; %scales up ground truth positions

R = factor*[zeros(2,1) rand(2,7);zeros(1,nMics)];
R(1,1)=rand(1,1);

for kk=1: (size(R,2)/2)
    phi=2*pi*rand;
    Rot=[cos(phi) sin(phi); -sin(phi) cos(phi)];
    R(:,kk*2)=R(:,kk*2-1) + [Rot*rigDist*[1; 0] ; 0] ;
end

%%Sanity check on receiver ditances = rigDist
norm(R(:,1)-R(:,2))
norm(R(:,3)-R(:,4))
norm(R(:,5)-R(:,6))
norm(R(:,7)-R(:,8))

S = factor*[0 0;rand(1,1)*0.2 1+rand(1,1);rand(1,2)];
d=toa_calc_d_from_xy(R,S);

if findIntegerMeasurements
   %find gt so that d is integer
   dRound=round(d);
   [I,J,D]=find(dRound);
   
   [R,S]=bundleToaMicRigsUncalibMics2DSounds3D(D,I,J,R(1:2,:),S,0.01)
   R=[R ; repmat(1,size(R,2))];
    
end

RpGt = R(2,:);  %Receiver Porjection Grund Truth
SpGt = S(2,:);  %Sender Porjection Grund Truth



d2 = d.^2;

d11 = d2(1,1);
d1j = d2(1,:)-d2(1,1);
di1 = d2(:,1)-d2(1,1);
dij = d2-repmat(d2(:,1),1,2)-repmat(d2(1,:),nMics,1)+d2(1,1);
Rtilde = [0 dij(2:end,2)'];
Stilde = [0 1];

l = Rtilde/RpGt;
b = mean(-2*SpGt/(l))-mean(Stilde) %So s_(2,1)=l*b/-2 and s_(2,2)=l*(Stilde(2)+b)/-2
RpGt
inv(l)*Rtilde
SpGt 
(l*(Stilde+b))/(-2)
h = 1/(l^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% so far so good

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
v7 = R(1,7)^2; 
v8 = R(1,8)^2;
u12 = R(1,1)*R(1,2);
u34 = R(1,3)*R(1,4);
u56 = R(1,5)*R(1,6);
u78 = R(1,7)*R(1,8);



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
%testa di(7)-
%di(8)-

clear A B

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
xxx = [h b v1 v2 v3 v4 v5 v6 v7 v8 u12 u34 u56 u78]'; %yeah, derp.


%%%%%%%%%%%%%%%%
%% ridDist^2 -equations
% We now use the three equations ||r_1-r_2||=||r3-r4||=||r5-r6||=||r7-r8||


% ||r_1-r_2||-||r3-r4||=0 foloows here
 norm(R(:,2)-R(:,1))
%U12^2
 %groud truth test, followed by Parametrization ground test. Each term
 %above and below should be equal.
R(:,2)'*R(:,2)      -2*R(:,2)'*R(:,1) + R(:,1)'*R(:,1) - R(:,3)'*R(:,3)      + 2*R(:,3)'*R(:,4)                 - R(:,4)'*R(:,4) %ok
v2 + Rtilde(2)^2*h  -2* u12           + v1             -  v3 - Rtilde(3)^2*h + 2*Rtilde(3)*Rtilde(4)*h + 2*u34  - v4-Rtilde(4)^2*h  %ok

B(8,1) = 0; %we take the two rig distances minus each other 
A(8,1:14)=[(Rtilde(2)^2-Rtilde(3)^2 + 2*Rtilde(3)*Rtilde(4) - Rtilde(4)^2) 0 1 1 -1 -1 0 0 0 0 -2 2 0 0];

%%%%%%%%%%%%%%%%
%% 

% ||r_1-r_2||-||r5-r6||=0 foloows here
 norm(R(:,6)-R(:,5)); %should be =rigDist

 %groud truth test, followed by Parametrization ground test. Each term
 %above and below should be equal.
R(:,2)'*R(:,2)      -2*R(:,2)'*R(:,1)   + R(:,1)'*R(:,1)  - R(:,5)'*R(:,5)       + 2*R(:,5)'*R(:,6)                  - R(:,6)'*R(:,6) %ok

v2 + Rtilde(2)^2*h  -2* u12             + v1              - v5 - Rtilde(5)^2*h   + 2*Rtilde(5)*Rtilde(6)*h + 2*u56   - v6-Rtilde(6)^2*h  %ok

B(9) = 0;
A(9,1:14)=[(Rtilde(2)^2-Rtilde(5)^2 + 2*Rtilde(5)*Rtilde(6) - Rtilde(6)^2) 0 1 1 0 0 -1 -1  0 0 -2 0 2 0];


%%%%%%%%%%%%%%%%
%% 

% ||r_1-r_2||-||r7-r8||=0 foloows here
 norm(R(:,8)-R(:,7)); %should be =rigDist

 %groud truth test, followed by Parametrization ground test. Each term
 %above and below should be equal.
R(:,2)'*R(:,2)      -2*R(:,2)'*R(:,1)   + R(:,1)'*R(:,1)  - R(:,7)'*R(:,7)       + 2*R(:,7)'*R(:,8)                  - R(:,8)'*R(:,8) %ok

v2 + Rtilde(2)^2*h  -2* u12             + v1              - v7 - Rtilde(7)^2*h   + 2*Rtilde(7)*Rtilde(8)*h + 2*u78   - v8-Rtilde(8)^2*h  

B(10) = 0;
A(10,1:14)=[(Rtilde(2)^2-Rtilde(7)^2 + 2*Rtilde(7)*Rtilde(8) - Rtilde(8)^2) 0 1 1 0 0 0 0 -1 -1  -2 0 0 2];


%% %%%%%%%%%%%%%%

[A*xxx B]

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

if 1
%     C1 = ceil(100*rand(size(C1)));
%     C1 = ceil(100*rand(size(C1)));
%     id1 = find((C1~=0) .* (C1 ~= 1));
%     id2 = find((C2~=0) .* (C2 ~= 1));
%     C1(id1) = ceil(100*rand(length(id1),1));
%     C2(id2) = ceil(100*rand(length(id2),1));
    C1 = round(100*C1);
    C2 = round(100*C2);
end
vu = C1*monvek+C2; %the other monomials expresad as al incomb of h, b,v u12
%Vu thus consists of the unknowns [v2 v3 v4 v5 v6 v7 v8 u34 u56 u78] inthat
%order, express in the unknowns h b v1 u12
%% THe equations to be solved by polysolve
% We here use the extra equations v1*v2=u12^2, v3*v4=u34^2, v5*v6=u56^2, v7*v8=u78^2 
eqs(1,1)=monvek(3)*vu(1)-monvek(4)^2;  %v1*v2=u12^2
eqs(2,1)=vu(2)*vu(3)-vu(8)^2; %v3*v4=u34^2
eqs(3,1)=vu(4)*vu(5)-vu(9)^2; %v5*v6=u56^2
eqs(4,1)=vu(6)*vu(7)-vu(10)^2; %v7*v8=u78^2 
gt = xxx(unknownsInd);
norm([evaluate(vu,gt) - xxx(extractedColsFromAInd)])/norm(xxx(extractedColsFromAInd))
evaluate(eqs,gt)

eqs2m2(100000*eqs,'eqs_8_2.m2')
%eval(['!m2<','eqs_8_2.m2'])
eqs2 = expandtodegree(eqs,[1 1 1 1],3);


for i = 1:length(eqs)
    [xx,yy] = polynomials2matrix(eqs(i));
    yy = monvec2matrix(yy);
    zz{i} = yy;
end


% eqs2 = generate_equations(eqs,3);

settings.dim = 6;
%settings.basis_size = 14;
settings.basis_size = 10;
settings.action_variable = 1;
[sols stats] = polysolve(eqs2, settings);

[err,errlist] = evaluate_solutions(sols,gt,settings);

err
