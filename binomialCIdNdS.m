function [u, l]= binomialCIdNdS(N, S, expected) 
    


t=N+S; %number of trials

e=N;
p=N/t;



[~,b] = binofit(N, t);

l=b(1)*t;
u=b(2)*t;

u=u/(t-u)/expected;
l=l/(t-l)/expected;



end
