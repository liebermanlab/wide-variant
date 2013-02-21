function [u, dNdS, l] = div_calc_dnds_ci(N,S, expected)

T=N+S;
p=N./T;
sd=sqrt(p.*(1-p)./T);
upper = p+(1.96*sd); 
lower = p-(1.96*sd); 

u=log(upper./(1-upper))/log(expected);
l=log(lower./(1-lower))/log(expected);

dNdS=log(N./S)./log(expected);


end
