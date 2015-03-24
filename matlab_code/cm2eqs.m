function eqs = cm2eqs(c,m);

eqs = multipol(0,zeros(size(m,1),1))*ones(size(c,1),1);
for kk = 1:size(c,1);
    cc = c(kk,:);
    jsel = find(cc);
    eqs(kk,1)=multipol(c(kk,jsel),m(:,jsel));
end;

