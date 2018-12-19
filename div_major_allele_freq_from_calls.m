function [MAF, majorNT,minorNT] = div_major_allele_freq_from_calls(ntcalls, ancnti)

%if second argument is provided, do MUTATION allele frequency

%same thing as div_major_allele_freq but takes in a 2d dimensional  matrix of calls
%(positions x isolates)

%built to take maNT as input-- could take anything

if ismember(ntcalls(1),'NATCGN')
    [~,ntcalls]=ismember(ntcalls,'NATCG');
    ntcalls=ntcalls-1;
elseif ismember(ntcalls(1),'natcg')
    [~,ntcalls]=ismember(ntcalls,'natcg');
    ntcalls=ntcalls-1;
elseif ~ismember(ntcalls(1),[0 1 2 3 4])
    error('Input to div_major_allele_freq_from_calls must be 0-4 or ATCGN or atcgn')
end

c=[sum(ntcalls==1,2) sum(ntcalls==2,2) sum(ntcalls==3,2) sum(ntcalls==4,2)];

if nargin < 2
    
    [sorted, sortedpositions] = sort(c,2);
    maxcount = sorted(:,end);
    minorcount = sorted(:,end-1);
    
    MAF=double(maxcount)./sum(c,2);
    minorAF=double(minorcount)./sum(c,2);
    
    majorNT = squeeze(sortedpositions(:,end));
    minorNT = squeeze(sortedpositions(:,end-1));
    
    MAF=squeeze(MAF);
    minorAF=squeeze(minorAF);
    
else
    
    if ismember(ancnti(1),'NATCGN')
        [~,ancnti]=ismember(ancnti,'NATCG');
        ancnti=ancnti-1;
    elseif ismember(ntcalls(1),'natcg')
        [~,ancnti]=ismember(ancnti,'natcg');
        ancnti=ancnti-1;
    elseif ~ismember(ancnti(1),[0 1 2 3 4])
        error('Ancestor input to div_major_allele_freq_from_calls must be 0-4 or ATCGN or atcgn')
    end    
    
    total=sum(c,2);
    c(sub2ind(size(c),1:numel(ancnti),ancnti'))=0;
    
    MAF=sum(c,2)./total;
    MAF(total<2)=0;
    
end




return

