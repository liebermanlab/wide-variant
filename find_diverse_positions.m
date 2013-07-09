function [p, numfields, coveragethresholds] = find_diverse_positions(params, SampleDirs, SampleNames, showscatters, parallel, jobsubmitoptions)

global TEMPORARYFOLDER; 

%p is a structure that contains the list of genomic positions meetings
%params

%coveragethresholds -- contains cdf cutoffs .01:.01:1.0.
%only counts positions where at least 1 read aligned
%purpose of this data structure is to allow downstream processes to
%remove positions with excess coverage in a way that accounts for
%variation in coverage between samples



if parallel==1
    if ~exist('find_diverse_positions_log')
        mkdir('find_diverse_positions_log')
    end
    
    %run on control -- not parallel
    fprintf(1,['\nAnalyzing isogenic control (' SampleNames{1} ')\n']) ;
    
    [p_control, coveragethresholds_control, MAF_control, numfields] = find_diverse_positions_single_sample(params, SampleDirs{1}, SampleNames{1}, [0], TEMPORARYFOLDER);
    
    
    %run others
    parallel_params={};
    for i=2:size(SampleNames)
        parallel_params{end+1}={params, SampleDirs{i}, SampleNames{i}, MAF_control, TEMPORARYFOLDER};        
    end
    run_parallel_matlab_commands('find_diverse_positions_single_sample', parallel_params, jobsubmitoptions, 1);
    
    
    %load files
    p=zeros(numel(p_control),1);
    coveragethresholds=zeros(numel(coveragethresholds_control),numel(SampleNames));
    p=p+p_control;
    coveragethresholds(:,1)=coveragethresholds_control;
    
    for i=2:size(SampleNames)
        %http://www.vsoch.com/2010/11/loading-dynamic-variables-in-a-static-workspace-in-matlab/
        diverse=load([TEMPORARYFOLDER '/diverse_' SampleNames{i} '.mat']);
        p=p+diverse.p_sample;
        coveragethresholds(:,i)=diverse.coveragethresholds_sample;
       delete([TEMPORARYFOLDER '/diverse_' SampleNames{i} '.mat'])
    end
    delete([TEMPORARYFOLDER '/diverse_' SampleNames{1} '.mat'])
    

    %return some information
    p=find(p>0);
        
    
else




   %Currently does not call ...single_sample because single_sample doesn't
   %allow interactivity



    %unpack parameters
    minorfreqthreshold=params.minorfreqthreshold;
    %maxreads_perstrand=params.maxreads_perstrand;
    minreads_perstrand=params.minreads_perstrand;
    minreads_perstrand_per_allele=params.minreads_perstrand_per_allele;
    min_bq=params.min_bq;
    min_mq=params.min_mq;
    max_td=params.max_td;
    min_td=params.min_td;
    max_sbp=params.max_sbp;
    max_bqp=params.max_bqp;
    %max_mqp=params.max_mqp;
    max_tdp=params.max_tdp;
    max_percent_indels=params.max_percent_indels;


    num_reqs=9;


    p=zeros(10^7,1);
    fprintf(1,'Finding diverse locations \n') ;


    coveragethresholds=zeros(100,size(SampleDirs));


    %Find diverse locations
    for i=1:length(SampleDirs)


        %load data
        fprintf(1,'Sample: %g  \n',i) ;
        load([SampleDirs{i} '/diversity.mat']);


        %coveragethresholds -- contains cdf cutoffs .01:.01:1.0.
        %only counts positions where at least 1 read aligned
        %purpose of this data structure is to allow downstream processes to
        %remove positions with excess coverage in a way that accounts for
        %variation in coverage between samples
        cov=sum(double(data(1:8,:)));
        cov(cov<1)=[]; cov=sort(cov);
        cutoffs=.01:.01:1;
        for j=1:numel(cutoffs);
            coveragethresholds(j,i)=cov(floor(cutoffs(j)*numel(cov)));
        end





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
        majorbqf=data(sub2ind(size(data),n1+8,positionsv))';
        majorbqr=data(sub2ind(size(data),n1+12,positionsv))';
        minorbqf=data(sub2ind(size(data),n2+8,positionsv))';
        minorbqr=data(sub2ind(size(data),n2+12,positionsv))';
        majormqf=data(sub2ind(size(data),n1+16,positionsv))';
        majormqr=data(sub2ind(size(data),n1+20,positionsv))';
        minormqf=data(sub2ind(size(data),n2+16,positionsv))';
        minormqr=data(sub2ind(size(data),n2+20,positionsv))';

        %Changed December 2012 to require forward and reverse strand qualities
        %   majorbq=(((f1.*data(sub2ind(size(data),n1+8,positionsv)))+(r1.*data(sub2ind(size(data),n1+12,positionsv))))./(f1+r1))';
        % majormq=(((f1.*data(sub2ind(size(data),n1+16,positionsv)))+(r1.*data(sub2ind(size(data),n1+20,positionsv))))./(f1+r1))';
        %  minorbq=(((f2.*data(sub2ind(size(data),n2+8,positionsv)))+(r2.*data(sub2ind(size(data),n2+12,positionsv))))./(f2+r2))';
        % minormq=(((f2.*data(sub2ind(size(data),n2+16,positionsv)))+(r2.*data(sub2ind(size(data),n2+20,positionsv))))./(f2+r2))';
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
        Treads= readsf > minreads_perstrand & readsr > minreads_perstrand & (f2' > minreads_perstrand_per_allele) & (r2' > minreads_perstrand_per_allele);
        % max reads per strand removed for now in favor of thresholding later
        % & readsf < maxreads_perstrand & readsr < maxreads_perstrand ;


        %Changed December 2012 to require forward and reverse strand qualities
        Tbq= ((majorbqf > min_bq) & (minorbqf > min_bq) & (majorbqr > min_bq) & (minorbqr > min_bq));
        Tmq = ((majormqf > min_mq) & (minormqf > min_mq) & (majormqr > min_mq) & (minormqr > min_mq));
        % Tbq= (majorbq > min_bq) & (minorbq > min_bq);
        % Tmq = (majormq > min_mq) & (minormq > min_mq);
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
        %fprintf(1,'Max requirements met: %g  \n',num_reqs) ;
        %num_reqs=max(allreqs);

        fprintf(1,'Positions meeting all requirements: %g  \n',sum(allreqs==num_reqs)) ;

        p(allreqs==num_reqs)=1;


        %Report how many positions are removed because of Isogenic control
        if i==1
            controlMAF=maf;
        else
            good=div_single_sample_test_thresholds(data, params, controlMAF, coveragethresholds(:,i));
            removedbycontrol=div_single_sample_test_thresholds(data, params, ones(size(controlMAF)), coveragethresholds(:,i));
            removedbycontrol(good>0)=0;
            fprintf(1,'Positions only removed from looseparameters because of Isogenic control: %g  \n',sum(removedbycontrol))
        end



        %Clickable scatter vs control, only plotting  10% of 'bad' points
        %  toplot=(rand(size(allreqs)) > .9);%plot one in every 10 points
        if showscatters > 0 & i~=1

            figure(130+i); clf; hold on;
            search_ind = find(removedbycontrol) ;
            title(SampleNames(i))
            plot(controlMAF,maf,'o','MarkerSize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); %, 'ButtonDownFcn',@clicked) ;
            plot(controlMAF(good>0),maf(good>0),'o','MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k'); %, 'ButtonDownFcn' ,@clicked) ;
            plot(controlMAF(removedbycontrol>0),maf(removedbycontrol>0),'o','MarkerSize', 5, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'ButtonDownFcn' ,@clicked) ;
            axis([.5 1 .5 1])

            disp('Paused. Hit any key continue...')
            pause
        end

    end


    %Find ll locations where at least one strain was diverse

    p=find(p);


    %Convert positions to Chr, Pos

    numfields=size(data,1);



end


        function clicked(src,event)

            %find data point clicked, this region written by Roy Kishony March
            %2012
            ac=get(gca,'position') ;
            fc=get(gcf,'position') ;
            pc=get(gcf,'CurrentPoint') ;
            xl=xlim ;
            yl=ylim ;
            ax = [fc(3) * [ac(1), ac(1)+ac(3)]]  ;
            ay = [fc(4) * [ac(2), ac(2)+ac(4)]]  ;
            x0 = xl(1) + (xl(2)-xl(1))*(pc(1)-ax(1))/(ax(2)-ax(1)) ;
            y0 = yl(1) + (yl(2)-yl(1))*(pc(2)-ay(1))/(ay(2)-ay(1)) ;

            search_ind = find(removedbycontrol) ;

            [~,ind] = min((controlMAF(search_ind)-x0).^2+(maf(search_ind)-y0).^2) ;
            ind = search_ind(ind) ;
            div_display_thresholds(data(:,ind), params, controlMAF(ind), coveragethresholds(:,i))


            figure(80);clf;

            %Plot number of calls
            subplot(5,1,1); hold on;
            title(['Number of calls....  log p(no strand bias) = ' num2str(SBp(ind))]);
            a=squeeze(data(1:8,ind));
            bar(reshape([a; nan(4,size(a,2))],4,[])','stacked')
            legend('A','T','C','G')
            ylabel('Number of reads')


            %Plot call quality
            subplot(5,1,2); hold on;
            title(['Average call quality.... log p = ' num2str(BQp(ind))]);
            a=squeeze(data(9:16,ind));
            bar(reshape([a; nan(4,size(a,2))],4,[])','grouped')
            axis([0.5 3.5 0 50])
            legend('Aq','Tq','Cq','Gq')
            ylabel('Average Base Quality')

            %Plot mapping quality
            subplot(5,1,3); hold on;
            title(['Average mapping quality.... log p = ' num2str(MQp(ind))]);
            a=squeeze(data(17:24,ind));
            bar(reshape([a; nan(4,size(a,2))],4,[])','grouped')
            axis([0.5 3.5 0 50])
            legend('Am','Tm','Cm','Gm')
            ylabel('Average Mapping Quality')


            %Plot tail distance forward
            subplot(5,1,4); hold on;
            title(['Average tail distance....  log p = ' num2str(TDFp(ind))]);
            a=squeeze(data(25:32,ind));
            bar(reshape([a; nan(4,size(a,2))],4,[])','grouped')
            axis([0.5 3.5 0 50])
            legend('Atd','Ttd','Ctd','Gtd')
            ylabel('Average Tail Distance')

            %Plot tail distance reverse
            subplot(5,1,5); hold on;
            title(['Average tail distance.... log p = ' num2str(TDRp(ind))]);
            a=squeeze(data(25:32,ind));
            bar(reshape([a; nan(4,size(a,2))],4,[])','grouped')
            axis([0.5 3.5 0 50])

            legend('Atd','Ttd','Ctd','Gtd')
            ylabel('Average Tail Distance')

        end
end






