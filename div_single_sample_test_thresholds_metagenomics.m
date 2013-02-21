function [allthresholdsmet] = div_single_sample_test_thresholds_metagenomics(d, p)

num_reqs=6;


NPositions=size(d,2);

%unpack parameters
minorfreqthreshold=p.minorfreqthreshold;
maxreads_perstrand=p.maxreads_perstrand;
minreads_perstrand=p.minreads_perstrand;
min_bq=p.min_bq;
min_mq=p.min_mq;
min_td=p.min_td;
max_td=p.max_td;
max_sbp=p.max_sbp;
max_bqp=p.max_bqp;
%max_mqp=p.max_mqp;
max_tdp=p.max_tdp;
max_percent_indels=p.max_percent_indels;

%unpack relevant stats
[~, majorNT, minorNT] = div_major_allele_freq(d);


allthresholdsmet=zeros(NPositions,1);
positionsv=(1:NPositions)';

n1=majorNT';
n2=minorNT';

minorfreq=(d(sub2ind(size(d),n2,positionsv))+d(sub2ind(size(d),n2+4,positionsv)))'./sum(d(1:8,:));
readsf=sum(d(1:4,:));
readsr=sum(d(5:8,:));
f1 = d(sub2ind(size(d),n1,positionsv)); %major allele counts on forward strand
r1 = d(sub2ind(size(d),n1+4,positionsv)); %major allele counts on forward strand
f2 = d(sub2ind(size(d),n2,positionsv)); %major allele counts on forward strand
r2 = d(sub2ind(size(d),n2+4,positionsv)); %major allele counts on forward strand
majorbq=((f1.*d(sub2ind(size(d),n1+8,positionsv)))+(r1.*d(sub2ind(size(d),n1+12,positionsv))))./(f1+r1);
majormq=((f1.*d(sub2ind(size(d),n1+16,positionsv)))+(r1.*d(sub2ind(size(d),n1+20,positionsv))))./(f1+r1);
minorbq=((f2.*d(sub2ind(size(d),n2+8,positionsv)))+(r2.*d(sub2ind(size(d),n2+12,positionsv))))./(f2+r2);
minormq=((f2.*d(sub2ind(size(d),n2+16,positionsv)))+(r2.*d(sub2ind(size(d),n2+20,positionsv))))./(f2+r2);
majortdf=(d(sub2ind(size(d),n1+24,positionsv)));
minortdf=(d(sub2ind(size(d),n2+24,positionsv)));
majortdr=(d(sub2ind(size(d),n1+28,positionsv)));
minortdr=(d(sub2ind(size(d),n2+28,positionsv)));
percent_indels=double(d(end,:)./double(sum(d(1:8,:))+d(end,:)));
SBp=d(end-6,:);
BQp=d(end-5,:);
MQp=d(end-4,:);
TDFp=d(end-3,:);
TDRp=d(end-2,:);


%Find true/false of meeting thresholds
Tminor = minorfreq > minorfreqthreshold;
Treads= (readsf > minreads_perstrand) & (readsr > minreads_perstrand) & (readsf < maxreads_perstrand) & (readsr < maxreads_perstrand) ;
Tbq= ((majorbq > min_bq) & (minorbq > min_bq))';
Tmq = ((majormq > min_mq) & (minormq > min_mq))';
Ttd = ((majortdf > min_td) & (majortdf < max_td) & (majortdr < max_td) & (majortdr > min_td)...
    & (minortdf > min_td) & (minortdf < max_td) & (minortdr > min_td) & (minortdr < max_td))';
Tid = percent_indels < max_percent_indels;
TSBp = SBp < max_sbp;
TBQp = BQp < max_bqp;
%TMQp = MQp < max_mqp;
TTDp = (TDFp < max_tdp) & (TDRp < max_tdp);




allreqs= Tminor + Treads + Tbq + Tmq + Ttd + Tid + TSBp + TBQp +TTDp; % TMQp 
allthresholdsmet((allreqs>=num_reqs))=1;



end


