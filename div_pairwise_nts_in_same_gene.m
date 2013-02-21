function [p, g] = div_pairwise_nts_in_same_gene(genes)

%June 2012 Tami Lieberman
%Takes a list of gene numbers, corresponding to genes associated with each
%diverse nt, and finds the percent overlap

nump=numel(genes);
pairmatrix=(repmat(genes,1,nump)- repmat(genes',numel(genes),1))==0 - eye(numel(genes)); %same gene number

p=sum(pairmatrix(:))/(nump^2-nump); %both numerator and denominator would be divided by 2, they cancel out

g=repmat(genes,1,nump);

g=unique(g(pairmatrix)); %which genes were found more than once

end
