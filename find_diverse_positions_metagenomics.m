function [p, Positions, numfields] = find_diverse_positions_metagenomics(fname, params, SampleDirs, SampleNames, ChrStarts)


%unpack parameters
minorfreqthreshold=params.minorfreqthreshold;
maxreads_perstrand=params.maxreads_perstrand;
minreads_perstrand=params.minreads_perstrand;
min_bq=params.min_bq;
min_mq=params.min_mq;
max_td=params.max_td;
min_td=params.min_td;
max_sbp=params.max_sbp;
max_bqp=params.max_bqp;
%max_mqp=params.max_mqp;
max_tdp=params.max_tdp;
max_percent_indels=params.max_percent_indels;





p=zeros(10^7,1);
fprintf(1,'Finding diverse locations \n') ;


%Find diverse locations
for i=1:length(SampleDirs)
    
    
    %load data
    fprintf(1,'Sample: %g  \n',i) ;
    load([SampleDirs{i} '/' fname]);
    
    %Parse data
    NPositions=size(data,2);
    [maf, majorNT, minorNT] = div_major_allele_freq(data);
    positionsv=(1:NPositions)';
    n1=majorNT';
    n2=minorNT';
    
    minorfreq=(double(data(sub2ind(size(data),n2,positionsv)))+double(data((sub2ind(size(data),n2+4,positionsv)))))'./sum(double(data(1:8,:)));
    readsf=sum(double(data(1:4,:)));
    readsr=sum(double(data(5:8,:)));
    f1 = data(sub2ind(size(data),n1,positionsv)); %major allele counts on forward strand
    r1 = data(sub2ind(size(data),n1+4,positionsv)); %major allele counts on forward strand
    f2 = data(sub2ind(size(data),n2,positionsv)); %major allele counts on forward strand
    r2 = data(sub2ind(size(data),n2+4,positionsv)); %major allele counts on forward strand
    majorbq=(((f1.*data(sub2ind(size(data),n1+8,positionsv)))+(r1.*data(sub2ind(size(data),n1+12,positionsv))))./(f1+r1))';
    majormq=(((f1.*data(sub2ind(size(data),n1+16,positionsv)))+(r1.*data(sub2ind(size(data),n1+20,positionsv))))./(f1+r1))';
    minorbq=(((f2.*data(sub2ind(size(data),n2+8,positionsv)))+(r2.*data(sub2ind(size(data),n2+12,positionsv))))./(f2+r2))';
    minormq=(((f2.*data(sub2ind(size(data),n2+16,positionsv)))+(r2.*data(sub2ind(size(data),n2+20,positionsv))))./(f2+r2))';
    majortdF=(data(sub2ind(size(data),n1+24,positionsv)))';
    minortdF=(data(sub2ind(size(data),n2+24,positionsv)))';
    majortdR=(data(sub2ind(size(data),n1+28,positionsv)))';
    minortdR=(data(sub2ind(size(data),n2+28,positionsv)))';
    percent_indels=double(data(end,:))./double(sum(data(1:8,:))+double(data(end,:)));
    SBp=data(end-6,:);
    BQp=data(end-5,:);
    MQp=data(end-4,:);
    TDFp=data(end-3,:);
    TDRp=data(end-2,:);
    
    
    
    
    
    %Find true/false of meeting thresholds
    Tminor = minorfreq > minorfreqthreshold;
    Treads= readsf > minreads_perstrand & readsr > minreads_perstrand & readsf < maxreads_perstrand & readsr < maxreads_perstrand ;
    Tbq= (majorbq > min_bq) & (minorbq > min_bq);
    Tmq = (majormq > min_mq) & (minormq > min_mq);
    Ttd = (majortdF > min_td) & (majortdF < max_td) & (majortdR < max_td) & (majortdR > min_td)...
        & (minortdF > min_td) & (minortdF < max_td) & (minortdR > min_td) & (minortdR < max_td);
    Tid = percent_indels < max_percent_indels;
    TSBp = SBp < max_sbp;
    TBQp = BQp < max_bqp;
    %TMQp = MQp < max_mqp;
    TTDp = (TDFp < max_tdp) & (TDRp < max_tdp);
    
    
    
    %Report how many positions met each requirement
    fprintf(1,'MinorAlleleFreq: %g  \n',sum(Tminor)) ;
    fprintf(1,'Cov: %g  \n',sum(Treads)) ;
    fprintf(1,'SBp: %g  \n',sum(TSBp)) ;
    fprintf(1,'BQp: %g  \n',sum(TBQp)) ;
    %fprintf(1,'MQp: %g  \n',sum(TMQp)) ;
    fprintf(1,'TDp: %g  \n',sum(TTDp)) ;
    fprintf(1,'maxIndels: %g  \n',sum(Tid)) ;
    fprintf(1,'minBQ: %g  \n',sum(Tbq)) ;
    fprintf(1,'minMQ: %g  \n',sum(Tmq)) ;
    fprintf(1,'acceptableTD: %g  \n',sum(Ttd)) ;
    
        
    %Records positions that met all requirements
    allreqs= Tminor + Treads + Tbq + Tmq + Ttd + Tid + TSBp + TBQp + TTDp; %TMQp 
    num_reqs=max(allreqs);
    nontest_reqs=Tminor + Treads + Tbq + Tmq + Ttd + Tid;
    fprintf(1,'Max requirements met: %g  \n',num_reqs) ;
    fprintf(1,'Positions meeting all requirements: %g  \n',sum(allreqs==num_reqs)) ;
    
    p(allreqs==num_reqs)=1;
    
    %plot statistics
   % div_data_hists(data, nontest_reqs, params, ['Whole_genome_hist_' SampleNames(i).Sample])
   % div_data_cum_hists(data, nontest_reqs, params, ['Whole_genome_cum_hist_' SampleNames(i).Sample])
    
    
end

%Find locations where at least one strain was diverse

p=find(p);


%Convert positions to Chr, Pos
Positions = p2chrpos(p, ChrStarts)';

numfields=size(data,1);



end
