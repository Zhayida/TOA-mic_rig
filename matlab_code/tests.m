clearvars

m = 2;
n = 4;


options.dim = 3;
options.origin = 1;

[data] = generate_mic_rig(m,n,options);

nedge = 18;
nnode = 8;

nmon  = nedge + nnode - 1 - (n+1);

neq   = nedge;


CC = zeros(neq,nmon);

rr = [];

for i = 1:m*2
    gg{i} = [i*ones(m*2,1) m*2+(1:n)'];
    rr    = [rr;gg{i}];
end

for i = 1:m
    rr = [rr;[2*i-1 2*i]];
end


rrr  = rr;

x    = [data.r,data.s];

C2gt = x(:,2:end)'*x(:,2:end);


k = 0;

for i = 1:length(rr);
    if (rr(i,1)==1)
        CC(i,rr(i,2)-1) = 1;
        mongt(rr(i,2)-1) = C2gt(rr(i,2)-1,rr(i,2)-1);
    else
        CC(i,rr(i,1)-1) = 1;
        CC(i,rr(i,2)-1) = 1;
        k = k+1;
        CC(i,nnode-1+k) = -2;
        mongt(nnode-1+k) = C2gt(rr(i,1)-1,rr(i,2)-1);

    end
end

for i = 1:nnode-1
    mongt(i) = C2gt(i,i);
end

d2 = data.d.^2;
d2 = d2';
d2 = d2(:);
d2 = [d2;data.l^2*ones(m,1)];

norm(CC*mongt'-d2)


%%
ee  = CC\d2;
bas = null(CC); 

k  = nmon - neq;
E  = eye(k);
for i = 1:k
    xv(i) = multipol(1,E(:,i));
end
xv = xv(:);
one = multipol(1,zeros(k,1));
zero = multipol(0,zeros(k,1));

monvec = ee*one + bas*xv;

xvgt = (mongt'-ee)'*bas;


rrcc = rrr;
rrcc = rrcc(rr(:,1)>1,:)-1;

CC2 = zero*zeros(nnode-1);
for i = 1:length(rrcc);
    kk = rrcc(i,:);
    CC2(kk(1),kk(2)) = monvec(nnode-1+i);
    CC2(kk(2),kk(1)) = monvec(nnode-1+i);
end

for i = 1:nnode-1
    CC2(i,i) = monvec(i);
end



%%
rr = nchoosek(1:nnode-1,4);
cc = nchoosek(1:nnode-1,4);

kkk = 1;
% eqs = multipol(length(rr)*length(cc));
for i = 1:length(rr);
    i
    for j = 1:length(cc)
        dd = detv(CC2(rr(i,:),cc(j,:)));
        if ~isempty(coeffs(dd))
        eqs(kkk) = dd;
        kkk = kkk+1;
        end
    end
end

%%
C = polynomials2matrix(eqs);
for i = 1:length(eqs);
    eqs(i) = eqs(i)/norm(C(i,:));
end

%%
settings.dim = 4;
settings.basis_size = 6;
settings.action_variable = 1;
eqs3  = generate_equations(eqs,1);
[sols,stats] = polysolve(eqs3,settings);
 