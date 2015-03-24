load settings_4_5
options.dimr = 3;
options.dims = 3;
m = 2;
n = 5;

options.dim = 3;
options.origin = 1;

[data] = generate_mic_rig(m,n,options);

[rpos,spos,flag] = solver_4_5(data);

%%
pts_gt = [data.r,data.s];
for i = 1:length(rpos);
    pts = [rpos{i},spos{i}];
    [err,pts_T,T] = transform(pts_gt,pts);
    err
end




% stats