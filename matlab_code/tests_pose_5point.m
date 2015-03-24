clearvars

m = 1;
n = 5;


options.dim = 3;
options.origin = 1;
options.rsyn = 1;
options.ssyn = 0;
[data] = generate_mic_rig(m,n,options);

f = data.f;
f2 = f.^2;
v = [data.s_offset;data.s_offset.^2;data.r(:,1)'*data.r(:,1);data.r(:,2)'*data.r(:,2);sum(data.r(:,1).*data.r(:,2));data.r(:,1);data.r(:,2);1];

M = zeros(2*m*n +  m, 2*n + 3 + 2*m*3 + 1);
k = 1;
for i = 1:2*m
    for j = 1:n
        M(k,j)   =  2*f(i,j);
        M(k,j+n) =  -1;
        M(k,2*n+i) = 1;
        M(k,end) = -(f2(i,j)) + data.s(:,j)'*data.s(:,j);
        M(k,(3+2*n+(i-1)*3+(1:3))) = -2*data.s(:,j);
        k  = k+1;
    end
end

for i = 1:m
    M(end,2*n+1) = 1;
    M(end,2*n+2) = 1;
    M(end,2*n+3) = -2;
    M(end,end)   = -data.l^2;
end

% normalization
for i = 1:size(M,1);
    M(i,:) = M(i,:)/norm(M(i,:));
end

% null 
bas = null(M(:,1:end-1));
e   = -M(:,1:end-1)\M(:,end);

gt  = (v(1:end-1)-e)'*bas;
gt  = gt';

% equations
E = eye(8);
for i = 1:length(E);
    xv(i) = multipol(1,E(:,i));
    
end

one = multipol(1,zeros(8,1));
zero = multipol(0,zeros(8,1));
;

if 0
e = round(randn(19,1)*20);
bas = round(randn(19,8)*20);
end
vv = e*one + (one*bas)*xv(:);

k = 1;
for i = 1:n
   eqs(k) = vv(i)*vv(i) - vv(n+i);
   k = k+1;
end

eqs(k) = vv(2*n+3+(1:3))'*vv(2*n+6+(1:3)) - vv(2*n+3);
eqs(k+1) = vv(2*n+3+(1:3))'*vv(2*n+3+(1:3)) - vv(2*n+1);



eqs(k+2) = vv(2*n+6+(1:3))'*vv(2*n+6+(1:3)) - vv(2*n+2);


eqs2m2(eqs,'pos_5point.m2');



%%
[cc,dd] = get_degree(eqs);

%%
for i = 1:length(eqs);
    [~,mmm{i}] = polynomials2matrix(eqs(i));
    mmm{i} = monvec2matrix(mmm{i});
end

nn = 4;
mm = [1,1,1,1,1,1,1,1];
multimons = monvec(create_upto_degree(mm,nn));

eqs3 = eqs(:)*multimons(:)';
eqs3 = eqs3(:);

[C,mon] = polynomials2matrix(eqs3);
size(C)

settings.dim = 40;
settings.basis_size = 100;
settings.action_variable = 1;
[sols,stats] = polysolve(eqs3,settings);