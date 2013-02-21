function callable = div_genomic_positions_with_potential_to_call_diversity(dirs, SampleNames, params,coveragethresholds, glength, parallel, submitoptions)


%Tami Lieberman December 2012
%Built based on find_diverse_positions.m



if parallel==1
    if ~exist('find_genomic_positions_with_potential_to_call_diversity_log')
        mkdir('find_genomic_positions_with_potential_to_call_diversity_log')
    end
    
    
    %run others
    parallel_params={};
    for i=1:size(SampleNames)
        parallel_params{end+1}={dirs{i}, SampleNames(i).Sample, params, coveragethresholds(:,i), glength};        
    end
    run_parallel_matlab_commands('div_genomic_positions_with_potential_to_call_diversity_single_s', parallel_params, submitoptions, 1);
    
    
    %load files
    callable=zeros(glength,numel(dirs));

    for i=1:size(SampleNames)
        %http://www.vsoch.com/2010/11/loading-dynamic-variables-in-a-static-workspace-in-matlab/
        c=load(['callable_' SampleNames(i).Sample '.mat']);
        callable(:,i)=c.callable_sample;
        delete(['callable_' SampleNames(i).Sample '.mat'])
    end

    
else

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


    callable=zeros(glength,size(dirs));
    fprintf(1,'Finding callable locations \n') 


    %Find diverse locations
    for i=1:length(dirs)

        
        maxreads_perstrand=coveragethresholds(params.maxreads_perstrand_percentile,i);
        %load data
        fprintf(1,'Sample: %g  \n',i) ;
        load([dirs{i} '/diversity.mat']);




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
        fprintf(1,'Cov: %g  \n',sum(Treads)) ;
        fprintf(1,'maxIndels: %g  \n',sum(Tid)) ;
        fprintf(1,'minBQ: %g  \n',sum(Tbq)) ;
        fprintf(1,'minMQ: %g  \n',sum(Tmq)) ;
        fprintf(1,'acceptableTD: %g  \n',sum(Ttd)) ;
        fprintf(1,'max_percentends: %g  \n',sum(Te)) ;


        %Records positions that met all requirements
        allreqs= Treads + Tbq + Tmq + Ttd + Tid + Te;
        fprintf(1,'Positions meeting all requirements: %g  \n',sum(allreqs==num_reqs)) ;

        callable(allreqs==num_reqs,i)=1;


    end

end



return

