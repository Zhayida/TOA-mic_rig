function [rpos,spos,flag] = solver_4_5(data,settings)

if nargin == 1;
    load settings_4_5
end

m = 2;
n = 5;

options.dim = 3;
options.origin = 1;

D = data.d;

Cm = [-ones(1,m*2-1);eye(m*2-1)];
Cn = [-ones(1,n-1);eye(n-1)];

DD = Cm'*D.^2*Cn;

[uu,ss,vv] = svd(DD);

R = (uu(:,1:3)*ss(:,1:3))';
S = vv(:,1:3)';

Rt = [zeros(3,1) R];
St = [zeros(3,1) S];
St = St/(-2);
Rt = Rt;



S = -S/2;

E = eye(11);
one = multipol(1,zeros(11,1));
zero = multipol(0,zeros(11,1));

D2 = D.^2;

for i = 1:11
    xv(i) = multipol(1,E(:,i));
end

H = [xv(1) xv(2) xv(3);
    xv(2) xv(4) xv(5);
    xv(3) xv(5) xv(6)];
b = [xv(7) xv(8) xv(9)]';

detH = detv(H);
adH = adjv(H);

l = xv(end-1);

if 0
    R = ceil(20*randn(size(R)));
    S = ceil(20*randn(size(S)));
    D2 = ceil(20*randn(size(D2)));
    data.l = 2;
end

eqs(1) = detv(H)*D2(1,1) - b'*adH*b;

k = 2;
for i = 1:n-1
    eqs(k) = detv(H)*(D2(1,i+1) - D2(1,1)) - S(:,i)'*adH*S(:,i) - 2*b'*adH*S(:,i);
    k      = k + 1;
end

for i = 1:2*m-1
    eqs(k) = (D2(i+1,1) - D2(1,1)) - R(:,i)'*H*R(:,i) + 2*b'*R(:,i);
    k      = k + 1;
end

eqs(k) = R(:,1)'*H*R(:,1) - l;

eqs(k+1) = R(:,2)'*H*R(:,2) - 2*R(:,2)'*H*R(:,3) + R(:,3)'*H*R(:,3) - l;

eqs(k+2) = detH*xv(11) - 1;


C = polynomials2matrix(eqs);



for i = 1:length(eqs)
    eqs(i) = eqs(i)/norm(C(i,:));
    
    [ccc,mmm{i}] = polynomials2matrix(eqs(i));
    mmm{i}       = monvec2matrix(mmm{i});
end


%%
[C,mon] = polynomials2matrix(eqs(end-5:end-1));
mm      = monvec2matrix(mon);

rrid = [1,2,7,8,9];
mmid = setdiff(1:10,rrid);


C = C(:,[rrid,[mmid,11]]);

mm = mm(:,[rrid,[mmid,11]]);


C = C(:,1:5)\C;
Mc = C(:,6:end);

eqsnew = matrix2polynomials(C(:,6:end),[mm(mmid,6:end);zeros(1,6)]);

for i = 1:5
    yv(rrid(i)) = -eqsnew{i};
end

E = eye(6);
for i = 1:length(mmid)
    yv(mmid(i)) = multipol(1,E(:,i));
end

Hn = [yv(1) yv(2) yv(3);
    yv(2) yv(4) yv(5);
    yv(3) yv(5) yv(6)];
bn = [yv(7); yv(8); yv(9)];

adHn = adjv(Hn);
detHn = detv(Hn);

eqsn(1) = multipol(1,E(:,6))*D2(1,1) - bn'*adHn*bn;

k = 2;
for i = 1:n-1
    eqsn(k) = multipol(1,E(:,6))*(D2(1,i+1) - D2(1,1)) - S(:,i)'*adHn*S(:,i) - 2*bn'*adHn*S(:,i);
    k      = k + 1;
end

eqsn(k) = detHn - multipol(1,E(:,6));

[mm,nn]  = get_degree(eqsn);


if (exist('settings','var') == 0)
NN = 7;
MM = 8;
m1max = [NN NN NN NN NN 1]';
tmax1 = MM;


for i = 1:length(eqsn(1:end));
    multimons{i} = monvec(create_upto_degree(m1max-mm{i},tmax1-nn(i)));
end
%
disp('gen template');
else
    
multimons = settings.multimons;
end
C = polynomials2matrix(eqsn);

eqs3 = [];
for i = 1:length(eqsn(1:end));
%     multimons{i} = monvec(create_upto_degree(m1max-mm{i},tmax1-nn(i)));
    eqsn(i) = eqsn(i)/norm(C(i,:));
    eqs3 = [eqs3(:);eqsn(i)*multimons{i}(:)];
end

[C1,mon1] = polynomials2matrix(eqs3);
mon1      = monvec2matrix(mon1);
id1       = find(mon1(6,:)==0);
id2       = find(mon1(6,:)==1);


%%
C2 = C1(:,[id1,id2]);
mon2 = mon1(:,[id1,id2]);

id1p = 1:length(id1);
id2p = (length(id1)+1):(length(id1)+length(id2));
[qq,rr,ee]=qr(C2(:,id1p));
zz =log10(abs(diag(rr/rr(1))))';

kk = min(find((zz)<-15));

%plot(zz)

eqs3 = cm2eqs(qq(:,(kk+1):end)'*C2(:,id2p),mon2(:,id2p));

mons3v = mon2-repmat([0 0 0 0 0 1]',1,size(mon2,2));
eqs3 = cm2eqs(qq(:,(kk+1):end)'*C2(:,id2p),mons3v(1:5,id2p));

if 1
    %%
    settings.basis_size = 50;
    settings.dim = 28;
    settings.action_variable = 1;
    
    %     eqs3 = generate_equations(eqsn,2);
    
    %     C = polynomials2matrix(eqs3);
    %     size(C)
    [sols,stats] = polysolve(eqs3,settings);

    
    sols_full = zeros(9,size(sols,2));
    sols_full(mmid,:) = sols;
    sols_full(rrid,:) = -Mc*[sols;ones(1,size(sols,2))];
    sols_full = sols_full(1:9,:);

    options.find_best = 0;
    [pos] = reconstruct_from_Cb(sols_full,Rt,St,3,options);
    
    
    if isempty(pos.x) == 0
        
        for i = 1:pos.n_good
            rpos{i} = pos.x{i};
            spos{i} = pos.y{i};
            flag = true;
        end
    else
        rpos = [];
        spos = [];
        flag = false;
    end
    
 
end
