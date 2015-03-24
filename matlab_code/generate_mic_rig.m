function [data] = generate_mic_rig(m,n,options)
% generate m microphone rigs and n speakers
% [data] = generate_mic_rig(m,n,options)

if ~isfield(options,'length')
    options.length = 0.2 + 0.1*rand;
end


if ~isfield(options,'dim')
    options.dim = 3;
end

if ~isfield(options,'origin')
    options.origin = 0;
end

if ~isfield(options,'rsyn'); % sync-rigs
    options.rsyn = 1;
end

if ~isfield(options,'ssyn'); % sync-senders
    options.ssyn = 1;
end

if ~isfield(options,'dimr'); % sync-senders
    options.dimr = options.dim;
end

if ~isfield(options,'dims'); % sync-senders
    options.dims = options.dim;
end

r1 = randn(options.dim,m);  % first end point of the rig
if options.origin
    r1(:,1) = zeros(options.dim,1); % set the first end point to be at the origin
end
nn = randn(options.dimr,m);  % direction of the rig
nn = nn./(ones(options.dimr,1)*sqrt(sum(nn.^2)));
ll = options.length;
r2 = r1 + ll*nn;            % second end point of the rig

r  = zeros(options.dimr,2*m);
r(:,1:2:2*m-1) = r1;
r(:,2:2:2*m)   = r2;

r = [r;zeros(3-options.dimr,2*m)];

s  = randn(options.dims,n);

s = [s;zeros(3-options.dims,n)];

data.d = compute_distance(r,s);

if options.ssyn == 0;
    s_offset = randn(1,n)';
    data.f   = data.d + (ones(2*m,1)*s_offset');
data.s_offset = s_offset;

end

if options.rsyn == 0;
    r_offset = randn(1,m)';
    data.f = data.d;
    for i = 1:m
        data.f(2*(i-1)+(1:2),:)   = data.d(2*(i-1)+(1:2),:) + (ones(2,n)*r_offset(i));
    end
    data.r_offset = r_offset;
end

data.l = ll;
data.s = s;
data.r1 = r1;
data.r2 = r2;
data.r  = r;



