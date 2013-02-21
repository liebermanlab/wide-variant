function p = strandbias(d)

[maf, m, n] = div_major_allele_freq(d);
%m is index (nt) of first major allele
%n is index (nt) of second major alle

%fishers exact test

f1 = d(sub2ind(size(d),m,1:numel(m))) ;
f2 = d(sub2ind(size(d),n,1:numel(n))) ;
r1 = d(sub2ind(size(d),m+4,1:numel(m))) ;
r2 = d(sub2ind(size(d),n+4,1:numel(n))) ;

p=zeros(1,numel(m));
for i=1:numel(m)
    r = fexact([f1(i) f2(i); r1(i) r2(i)], 0) ;
    p(i) = r(3);
end
