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
    div_data_hists(data, nontest_reqs, params, ['Whole_genome_hist_' SampleNames(i).Sample])
    div_data_cum_hists(data, nontest_reqs, params, ['Whole_genome_cum_hist_' SampleNames(i).Sample])
    
    
    
    %Clickable scatter vs control, only plotting  10% of 'bad' points
    good=(allreqs==num_reqs);
    toplot=(rand(size(allreqs)) > .9);%plot one in every 10 points
    if i==1
        controlMAF=maf;
    else
        figure(130+i); clf; hold on;
        search_ind = find(good) ;
        title(SampleNames(i).Sample)
        plot(controlMAF(~good & toplot),maf(~good & toplot),'o','MarkerSize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'ButtonDownFcn',@clicked) ;
        plot(controlMAF(good),maf(good),'o','MarkerSize', 5, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'ButtonDownFcn',@clicked) ;
    end
    
    
    
end

%Find ll locations where at least one strain was diverse

p=find(p);


%Convert positions to Chr, Pos
Positions = p2chrpos(p, ChrStarts)';

numfields=size(data,1);






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
        
        [~,ind] = min((controlMAF(search_ind)-x0).^2+(maf(search_ind)-y0).^2) ;
        ind = search_ind(ind) ;
        
        
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
