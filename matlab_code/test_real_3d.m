
clearvars
addpath ../../visionary
addpath ../../multipol

load './data/DMatchedByHand_G_Dim(mic,rec)=3,3_(nRigswith2mics,nSounds)=2,5_rowsAreMics_Mic1-2IsRig1_Mic3-4IsRig2_RigDistIs33.4cm.mat'

load './data/VisionReconstructionResultsVisionary_ G_Dim(mic,rec)=3,3_FIsrtFivePointsAreSoundSourcesSameOrderAsInD_LastFourPointsAreMicsSameOrderAsInD.mat'

pts = getpoints(s);
pts = pflat(pts);

%% RUnning the full 3d, 4 mics, 4 sounds calibrated rig case
rgt = pts(1:3,end-3:end)
soundsUsed=[1 2 3 5]  %4.1 5.6
sgt = pts(1:3,soundsUsed)

pts_gt = [rgt,sgt]


l = 0.334;
%%
m = 2;
n = 4;

D = D(:,soundsUsed);

Cm = [-ones(1,m*2-1);eye(m*2-1)];
Cn = [-ones(1,n-1);eye(n-1)];

DD = Cm'*D.^2*Cn;

[uu,ss,vv] = svd(DD);


dim = 3;




if 0
    
    R = (uu*ss)';
    S = vv';
    
    Rt = [zeros(3,1) R];
    St = [zeros(3,1) S];
    St = St/(-2);
    Rt = Rt;
    
else
    cutoff = 30;
    if (ss(1,1)/ss(dim,dim) < cutoff)
        xr11 = (uu(:,1:dim)*ss(1:dim,1:dim))';
        yr1 = vv(:,1:dim)';
    elseif ((ss(1,1)/ss(dim,dim) > cutoff && dim~=2))
        xr11 = (uu(:,1:dim))';
        yr1 = ss(1:dim,1:dim)*vv(:,1:dim)';
        
    elseif (dim==2)
        % This improve the numerics for 2D cases
        xr11 = (uu(:,1:dim)*sqrt(ss(1:dim,1:dim)))';
        yr1 = sqrt(ss(1:dim,1:dim))*vv(:,1:dim)';
    end
    
    if (ss(1,1)/ss(dim,dim) > cutoff)
    xr11 = (uu(:,1:dim)*sqrt(ss(1:dim,1:dim)))';
    yr1 = sqrt(ss(1:dim,1:dim))*vv(:,1:dim)';
    disp('type 1')
    else
    xr11 = DD';
    yr1  = eye(3);
    disp('type 2')
    end
    
    R = xr11;
    S = yr1;
    xtp = [zeros(dim,1) xr11];
    yt =[zeros(dim,1) yr1]/(-2);
    xt = xtp;
    
    St = yt;
    Rt = xtp;
    
    
    
end


S = -S/2;

E = eye(10);
one = multipol(1,zeros(10,1));
zero = multipol(0,zeros(10,1));

D2 = D.^2;

for i = 1:10
    xv(i) = multipol(1,E(:,i));
end

H = [xv(1) xv(2) xv(3);
    xv(2) xv(4) xv(5);
    xv(3) xv(5) xv(6)];
b = [xv(7) xv(8) xv(9)]';

detH = detv(H);
adH = adjv(H);

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

eqs(k) = R(:,1)'*H*R(:,1) - l^2;

eqs(k+1) = R(:,2)'*H*R(:,2) - 2*R(:,2)'*H*R(:,3) + R(:,3)'*H*R(:,3) - l^2;

eqs(k+2) = detH*xv(10) - 1;

eqs2m2(eqs(1:end),'eqsorg.m2');


%%
C = polynomials2matrix(eqs);



for i = 1:length(eqs)
    eqs(i) = eqs(i)/norm(C(i,:));
    
    [ccc,mmm{i}] = polynomials2matrix(eqs(i));
    mmm{i}       = monvec2matrix(mmm{i});
end



%%
[C,mon] = polynomials2matrix(eqs(end-5:end-1));
mm      = monvec2matrix(mon);

rrid = [1:2,7,8,9];
mmid = setdiff(1:9,rrid);



C = C(:,[rrid,[mmid,10]]);
mm = mm(:,[rrid,[mmid,10]]);


C = C(:,1:5)\C;

Mc = C(:,6:end);

% C(:,6:end) = ceil(randn(size(C(:,6:end)))*20);
% C(5,[6,7]) = 0;

% C(1:5,1:5) = eye(5);

eqsnew = matrix2polynomials(C(:,6:end),[mm(mmid,6:end);zeros(1,5)]);

for i = 1:5
    yv(rrid(i)) = -eqsnew{i};
end

E = eye(5);
for i = 1:length(mmid)
    yv(mmid(i)) = multipol(1,E(:,i));
end

Hn = [yv(1) yv(2) yv(3);
    yv(2) yv(4) yv(5);
    yv(3) yv(5) yv(6)];
bn = [yv(7); yv(8); yv(9)];

adHn = adjv(Hn);
detHn = detv(Hn);

eqsn(1) = multipol(1,E(:,5))*D2(1,1) - bn'*adHn*bn;

k = 2;
for i = 1:n-1
    eqsn(k) = multipol(1,E(:,5))*(D2(1,i+1) - D2(1,1)) - S(:,i)'*adHn*S(:,i) - 2*bn'*adHn*S(:,i);
    k      = k + 1;
end

eqsn(k) = detHn - multipol(1,E(:,5));

[mm,nn]  = get_degree(eqsn);

NN = 7;
MM = 8;
m1max = [NN NN NN NN 1]';
tmax1 = MM;
m2max = [NN NN NN NN 1]';
tmax2 = MM;

C = polynomials2matrix(eqsn);
eqs3 = [];
for i = 1:length(eqsn(1:end));
    multimons{i} = monvec(create_upto_degree(m1max-mm{i},tmax1-nn(i)));
    eqsn(i) = eqsn(i)/norm(C(i,:));
    eqs3 = [eqs3(:);eqsn(i)*multimons{i}(:)];
end

[C1,mon1] = polynomials2matrix(eqs3);
mon1      = monvec2matrix(mon1);
id1       = find(mon1(5,:)==0);
id2       = find(mon1(5,:)==1);


size(C1)
%%
C2 = C1(:,[id1,id2]);
mon2 = mon1(:,[id1,id2]);

id1p = 1:length(id1);
id2p = (length(id1)+1):(length(id1)+length(id2));
[qq,rr,ee]=qr(C2(:,id1p));
zz =log10(abs(diag(rr/rr(1))))';

kk = min(find((zz)<-15));

plot(zz)

eqs3 = cm2eqs(qq(:,(kk+1):end)'*C2(:,id2p),mon2(:,id2p));

mons3v = mon2-repmat([0 0 0 0 1]',1,size(mon2,2));
eqs3 = cm2eqs(qq(:,(kk+1):end)'*C2(:,id2p),mons3v(1:4,id2p));

%%
if 0
    
    %%
    C = polynomials2matrix(eqsn);
    
    eqs2m2(100000*eqsn,'eqs.m2');
    
    
    for i = 1:length(eqsn);
        eqsn(i) = eqsn(i)/norm(C(i,:));
    end
    
    gtn = [gt(mmid);det(Hgt)];
    
    evaluate(eqsn,gtn)
    
    C = polynomials2matrix(eqsn);
    for i = 1:length(eqsn)
        eqsn(i) = eqsn(i)/norm(C(i,:));
        
        [ccc,mmmm{i}] = polynomials2matrix(eqsn(i));
        mmmm{i}       = monvec2matrix(mmmm{i});
    end
    
    
    
    %% saturation
    [C,mon] = polynomials2matrix(eqsn);
    mon     = monvec2matrix(mon);
    id      = find(mon(5,:)==1);
    
end

if 1
    %%
    settings.basis_size = 20;
    settings.dim = 12;
    settings.action_variable = 1;

    [sols,stats] = polysolve(eqs3,settings);
    
    stats
    
    %%
    
    sols_full = zeros(9,size(sols,2));
    sols_full(mmid,:) = sols;
    sols_full(rrid,:) = -Mc*[sols;ones(1,size(sols,2))];
    
    
    options.find_best = 0;
    
    pos = reconstruct_from_Cb(sols_full,Rt,St,3,options);
    
    for i = 1:pos.n_good
        ptss{i} = [pos.x{i},pos.y{i}];
        [err,pts_T{i},T] = transform(ptss{i},pts_gt); %changed so that vision points alings t reconstructed, not other way around.
        errlist(i) = err
    end
    
    
    errlist
    [nouse,mid] = min(errlist)
    
    rrr = ptss{mid}(:,1:4);
    sss = ptss{mid}(:,5:end);
    %%
    good = pts_T{mid}
    ptss{mid}
    %
    figure(10);clf;hold on;
    rgt = good(:,1:4); %take out vision ground truth that is aligned with polysolve solution
    sgt = good(:,5:end);
    plot3(rrr(1,:),rrr(2,:),rrr(3,:),'rs','MarkerSize',12,'LineWidth',2); hold on;
    plot3(sss(1,:),sss(2,:),sss(3,:),'ro','MarkerSize',12,'LineWidth',2); hold on;
    

    plot3(rgt(1,:),rgt(2,:),rgt(3,:),'bs','MarkerSize',12,'LineWidth',2)
    plot3(sgt(1,:),sgt(2,:),sgt(3,:),'bo','MarkerSize',12,'LineWidth',2)
    legend('mic. - sound','speaker - sound','mic. - vision','speaker - vision')
    axis equal
     axis([-0.7 0.5 -0.6 0.5 -0.1 0.8])
    set(gca,'fontsize',15)
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    
        %plotting the rigs
    plot3(rrr(1,[1 2]),rrr(2,[1 2]),rrr(3,[1 2]),'r','MarkerSize',22,'LineWidth',1.5);
    plot3(rrr(1,[3 4]),rrr(2,[3 4]),rrr(3,[3 4]),'r','MarkerSize',22,'LineWidth',1.5);
    
    
    %set(gcf, {'markers'},{10})
    for i = 1:4
        for j = 1:4
            
            FF(i,j) = norm(rrr(:,i) - sss(:,j));
            
        end
    end

    FF - D
   
    RMSEMics        = sqrt(mean(sum((rrr-rgt).^2))) %Root Mean square distance error
    RMSESpreakers   = sqrt(mean(sum((sss-sgt).^2))) 
    %MEMics        = mean(sqrt(sum((rrr-rgt).^2))) %Mean distance error
    %MESpreakers   = mean(sqrt(sum((sss-sgt).^2))) 
    
    
end