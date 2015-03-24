function [err,pts_T,T] = transform(pts_gt,pts)

m = size(pts_gt,2);

[nouse,pts_T,T]=procrustes(pts_gt',pts');
pts_T = pts_T';
err=sqrt(sum(sum((pts_T-pts_gt).^2))/(m-1));

end