function [MAF, minorAF, majorNT, minorNT] = strains_major_allele_freq(calls)

nts='ATCG';
c=zeros(4,size(calls,1));
for p=1:size(calls,1)
    l=calls(p,:);
    for nt=1:4
        c(nt,p)=sum(l==nts(nt));
    end
end

%below based on div_major_allele_freq
[sorted, sortedpositions] = sort(c,1);
maxcount = sorted(end,:,:);
MAF=maxcount./sum(c,1);

majorNT = squeeze(sortedpositions(end,:,:));
minorNT = squeeze(sortedpositions(end-1,:,:));


minorcount = sorted(end-1,:,:);

minorAF=double(minorcount)./sum(c,1);
minorAF=squeeze(minorAF);

MAF=squeeze(MAF);

return

