function G = expandtodegree(G,Gdeg,maxdeg);
%maxdeg = 6;
gout = [];
for kk=1:length(G);
    gout = [gout;generate_equations(G(kk),maxdeg-Gdeg(kk))];
end
G = gout;
