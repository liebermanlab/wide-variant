function [u, l]= binomialCIs(s, n) 
    
% s is number of successes
% n is number of trials


e=s;
p=s/n;


sd=sqrt(n*p*(1-p));  
l=e-1.96*sd;   
u=e+1.96*sd;   

end
