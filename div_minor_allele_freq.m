function [minorAF] = div_minor_allele_freq(cnts)

c=cnts(1:4,:,:)+cnts(5:8,:,:);

[sorted, ~] = sort(c,1);
minorcount = sorted(end-1,:,:);

minorAF=double(minorcount)./sum(c,1);

minorAF=squeeze(minorAF);
return

