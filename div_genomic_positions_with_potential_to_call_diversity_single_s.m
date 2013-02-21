function callable = div_genomic_positions_with_potential_to_call_diversity_single_s(dir, samplename, params,coveragethresholds, glength, generatehists)


%unpack parameters
minreads_perstrand=params.minreads_perstrand;
min_bq=params.min_bq;
min_mq=params.min_mq;
max_td=params.max_td;
min_td=params.min_td;
max_percent_indels=params.max_percent_indels;
max_percent_ends=params.max_percent_ends;

%thresholdsunused
%minorfreqthreshold=params.minorfreqthreshold;
%max_sbp=params.max_sbp;
%max_bqp=params.max_bqp;
%max_mqp=params.max_mqp;
%max_tdp=params.max_tdp;


num_reqs=6;


callable=zeros(glength,1);

f=fopen(['find_genomic_positions_with_potential_to_call_diversity_log/' samplename '.txt'] ,'w');


maxreads_perstrand=coveragethresholds(params.maxreads_perstrand_percentile);
%load data

load([dir '/diversity.mat']);




%Parse data

readsf=sum(double(data(1:4,:)));
readsr=sum(double(data(5:8,:)));
bqf=(sum(double(data(9:12,:)).*double(data(1:4,:)))./readsf);
bqr=(sum(double(data(13:16,:)).*double(data(5:8,:)))./readsr);
mqf=(sum(double(data(17:20,:)).*double(data(1:4,:)))./readsf);
mqr=(sum(double(data(21:24,:)).*double(data(5:8,:)))./readsr);
tdF=(sum(double(data(25:28,:)).*double(data(1:4,:)))./readsf);
tdR=(sum(double(data(29:32,:)).*double(data(5:8,:)))./readsr);
percent_indels=double(data(end,:))./double(sum(data(1:8,:))+double(data(end,:)));
percent_ends=double(data(end-1,:))./sum(data(1:8,:));



%Find true/false of meeting thresholds
Treads= readsf > minreads_perstrand & readsr > minreads_perstrand & readsf < maxreads_perstrand & readsr < maxreads_perstrand ;

Tbq= ((bqf > min_bq) & (bqr > min_bq));
Tmq = ((mqf > min_mq) & (mqr > min_mq));
Ttd = ((tdF > min_td) & (tdF < max_td) & (tdR < max_td) & (tdR > min_td));
Tid = percent_indels < max_percent_indels;
Te = percent_ends < max_percent_ends;



%Report how many positions met each requirement
fprintf(f,'Cov: %g  \n',sum(Treads)) ;
fprintf(f,'maxIndels: %g  \n',sum(Tid)) ;
fprintf(f,'minBQ: %g  \n',sum(Tbq)) ;
fprintf(f,'minMQ: %g  \n',sum(Tmq)) ;
fprintf(f,'acceptableTD: %g  \n',sum(Ttd)) ;
fprintf(f,'max_percentends: %g  \n',sum(Te)) ;


%Records positions that met all requirements
allreqs= Treads + Tbq + Tmq + Ttd + Tid + Te;
fprintf(f,'Positions meeting all requirements: %g  \n',sum(allreqs==num_reqs)) ;

callable(allreqs==num_reqs)=1;

callable_sample=callable;


% plot statistics
if generatehists > 0
    
    %Get Tminor
    NPositions=size(data,2);
    [~, majorNT, minorNT] = div_major_allele_freq(data);
    positionsv=(1:NPositions)';
    n1=majorNT';
    n2=minorNT';
    minorfreq=(double(data(sub2ind(size(data),n2,positionsv)))+double(data((sub2ind(size(data),n2+4,positionsv)))))'./sum(double(data(1:8,:)));
    Tminor = minorfreq > minorfreqthreshold;
    nontest_reqs=Tminor + Treads + Tbq + Tmq + Ttd + Tid + Te;
    
    %Plot
    div_data_hists(data, nontest_reqs, params, ['Whole_genome_hist_' samplename])
end



save(['callable_' samplename], 'callable_sample')

end