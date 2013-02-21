function div_display_thresholds(d, p, controlMAF, coveragethresholds)

%This was written to perform on a giant matrix and adapted to report
%information for a single genomic position and strand by simply displaying 


%p here is the set of parameters

num_reqs=11;


NPositions=size(d,2);

%unpack parameters
minorfreqthreshold=p.minorfreqthreshold;
maxreads_perstrand=coveragethresholds(p.maxreads_perstrand_percentile);
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
max_percent_ends=p.max_percent_ends;
min_control_MAF=p.min_control_MAF;
minreads_perstrand_per_allele=p.minreads_perstrand_per_allele;

%unpack relevant stats
[~, majorNT, minorNT] = div_major_allele_freq(d);


allthresholdsmet=zeros(NPositions,1);
positionsv=(1:NPositions)';

n1=majorNT';
n2=minorNT';

minorfreq=double((d(sub2ind(size(d),n2,positionsv))+d(sub2ind(size(d),n2+4,positionsv))))'./sum(d(1:8,:));
readsf=sum(d(1:4,:));
readsr=sum(d(5:8,:));
f1 = d(sub2ind(size(d),n1,positionsv)); %major allele counts on forward strand
r1 = d(sub2ind(size(d),n1+4,positionsv)); %major allele counts on forward strand
f2 = d(sub2ind(size(d),n2,positionsv)); %minor allele counts on forward strand
r2 = d(sub2ind(size(d),n2+4,positionsv)); %minor allele counts on forward strand
majorbqf=d(sub2ind(size(d),n1+8,positionsv));
majorbqr=d(sub2ind(size(d),n1+12,positionsv));
minorbqf=d(sub2ind(size(d),n2+8,positionsv));
minorbqr=d(sub2ind(size(d),n2+12,positionsv));
majormqf=d(sub2ind(size(d),n1+16,positionsv));
majormqr=d(sub2ind(size(d),n1+20,positionsv));
minormqf=d(sub2ind(size(d),n2+16,positionsv));
minormqr=d(sub2ind(size(d),n2+20,positionsv));
majortdf=(d(sub2ind(size(d),n1+24,positionsv)));
minortdf=(d(sub2ind(size(d),n2+24,positionsv)));
majortdr=(d(sub2ind(size(d),n1+28,positionsv)));
minortdr=(d(sub2ind(size(d),n2+28,positionsv)));
percent_indels=double(d(end,:))./(sum(d(1:8,:))+double(d(end,:)));
percent_ends=double(d(end-1,:))./sum(d(1:8,:));

SBp=d(end-6,:);
BQp=d(end-5,:);
MQp=d(end-4,:);
TDFp=d(end-3,:);
TDRp=d(end-2,:);


%Find true/false of meeting thresholds
Tminor = minorfreq > minorfreqthreshold;
Treads= (readsf > minreads_perstrand) & (readsr > minreads_perstrand) & (readsf < maxreads_perstrand) & (readsr < maxreads_perstrand) & (f2' > minreads_perstrand_per_allele) & (r2' > minreads_perstrand_per_allele);
Tbq= ((majorbqf > min_bq) & (minorbqf > min_bq) & (majorbqr > min_bq) & (minorbqr > min_bq))';
Tmq = ((majormqf > min_mq) & (minormqf > min_mq) & (majormqr > min_mq) & (minormqr > min_mq))';
Ttd = ((majortdf > min_td) & (majortdf < max_td) & (majortdr < max_td) & (majortdr > min_td)...
    & (minortdf > min_td) & (minortdf < max_td) & (minortdr > min_td) & (minortdr < max_td))';
Tid = percent_indels < max_percent_indels;
Te = percent_ends < max_percent_ends;


TSBp = SBp < max_sbp;
TBQp = BQp < max_bqp;
%TMQp = MQp < max_mqp;
TTDp = (TDFp < max_tdp) & (TDRp < max_tdp);

if numel(controlMAF > 1)
    control = controlMAF > min_control_MAF;
else
    control = ones(Npositions,1);
end


allreqs= Tminor + Treads + Tbq + Tmq + Ttd + Tid + Te + TSBp + TBQp +TTDp+ control; % TMQp 
allthresholdsmet((allreqs==num_reqs))=1;

outcome={['yes'],['no']};

fprintf(1,['Frequency OK?  ' num2str(Tminor) '\n'])
fprintf(1,['Number of reads per strand in acceptable region?  ' num2str(Treads) '\n'])
fprintf(1,['Base quality OK?  ' num2str(Tbq) '\n'])
fprintf(1,['Mapping Quality OK?  ' num2str(Tmq) '\n'])
fprintf(1,['Tail Distance OK?  ' num2str(Ttd) '\n'])
fprintf(1,['Sufficiently low indel support in region?  ' num2str(Tid) '\n'])
fprintf(1,['Sufficiently low # of reads on ends?  ' num2str(Te) '\n'])
fprintf(1,['No strand bias?  ' num2str(TSBp) '\n'])
fprintf(1,['No base quality bias?  ' num2str(TBQp) '\n'])
fprintf(1,['No tail distance bias?  ' num2str(TTDp) '\n'])
fprintf(1,['Isogenic in isogenic control?  ' num2str(control) '\n'])




end


