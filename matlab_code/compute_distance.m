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