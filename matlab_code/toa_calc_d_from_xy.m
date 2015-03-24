function [ D ] = toa_calc_d_from_xy( rec, trans )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


m=size(rec,2);
n=size(trans,2);
D=zeros(m,n);

for kk=1:m 
    for jj=1:n
        D(kk,jj)=norm(rec(:,kk)-trans(:,jj));
    end
end


end

