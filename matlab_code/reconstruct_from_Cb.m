function [sols] = reconstruct_from_Cb(sols_full,xt,yt,d,options);

if nargin < 5;
    options.find_best = 1;
end
if isfield(options,'find_best')==0;options.find_best = 1;end
if isfield(options,'norm')==0;options.norm = 1;end

find_best = options.find_best;

kk = 1;
ok = 0;

if iscell(sols_full) && ~isempty(sols_full)
    nn = length(sols_full);
    sols_tmp = sols_full;
    sols_full = zeros(size(sols_tmp{1},1),nn);
    for i = 1:nn
        sols_full(:,i) = sols_tmp{i};
    end
end

if size(sols_full,1) > 5
    nc = 6;D = 3;
else
    nc = 3;D = 2;
end




for i = 1:size(sols_full,2);
    c = (sols_full(1:nc,i));
    if D == 3;
        CC{i} = [c(1) c(2) c(3) ; c(2) c(4) c(5) ; c(3) c(5) c(6)];
    else
        CC{i} = [c(1) c(2) ; c(2) c(3)];
    end
    try
        invC = inv(CC{i});
        l    = chol(invC);
        ok   = 1;
    catch % if not positive-definite, throw away
        %         disp(' C is not positive definite')
        l    = 0;
        continue;
    end
    C{kk} = CC{i};
    b{kk} = sols_full(end-D+1:end,i);
    L{kk} = l;
    solsvec{kk} = sols_full(:,i);
    kk = kk+1;
end
n_good = kk -1;
clear kk;

if ok
    %% reconstruct x, y
    err =  inf;
    % clear err;
    for kk = 1:length(L);

        x{kk} = inv(L{kk}')*(xt)*options.norm;
        y{kk} = L{kk}*bsxfun(@plus,yt,b{kk})*options.norm;

        if find_best
            dd = compute_distance(x{kk},y{kk});
            err(kk) = norm(d*options.norm - dd,'fro');
        end
    end
    
%     err
    % choose the one with smallest residue
    if find_best
        bestid = find(err == min(err));
        xbest   = x{bestid};
        ybest   = y{bestid};
        solsvec_best  = solsvec{bestid};
        errbest = min(err);
        sols.bestid = bestid;
    end
    sols.C = C;
    sols.L = L;
    sols.b = b;
    sols.x = x;
    sols.y = y;
    sols.solvec = solsvec;
    sols.n_good = n_good;
%     stats.n_real = n_real;
%     stats.n_good = length(L);
else
    sols = [];
    sols.n_good = 0;
    sols.solvec = [];
    sols.C = [];
    sols.L = [];
    sols.b = [];
    sols.x = [];
    sols.y = [];
%     disp('There is no solution found !');
    %     stats = [];
end

end


function [d] = compute_distance(x,y);
m = size(x,2);
n = size(y,2);
D = size(x,1);
% repmat [x(:,1),x(:,2),...; x(:,1),x(:,2)...;
vv = ones(n,m);
vv = ([1:m])'*ones(1,n);
vv = vv(:);
ind = (vv')+((1:length(vv))-1)*m;
ff  = zeros(m,m*n);
ff(ind) = 1;
xx = x*sparse(ff);

% repmat [y(:,1)...; y(:,2),...]
vv = ones(m,n);
vv = vv*(diag([1:n]));
vv = vv(:);
ind = vv'+((1:length(vv))-1)*n;
ff  = zeros(n,m*n);
ff(ind) = 1;
yy = y*sparse(ff);
d = sqrt( sum((xx-yy).^2));
d =reshape(d,m,n);
end