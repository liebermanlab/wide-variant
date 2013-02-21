
function div_clickable_scatter_sigcolor(x,y, xname, yname, sample, cnts, dwindows, iwindows, pos, annotations, RefGenome, ScafNames, SampleNames)


IsGenomeLoaded = false ;
scrsz = get(0,'ScreenSize');
strand=['g-','k-'];
window_size=floor(size(dwindows,1)/2);

figure;
clf ; hold on

sig=.01;
clrsig='rk';
clrcov='ck';

[maf, maNT, minorNT] = div_major_allele_freq(cnts);

%fishers exact test

sb=cnts(end-5,:,sample);
bq=cnts(end-4,:,sample);
mq=cnts(end-3,:,sample);
td=cnts(end-2,:,sample);
cov=sum(cnts(1:8,:,sample),1);

pt=.0001;
psig=(sb<pt | mq<pt | td<pt); %| bq<pt 
covsig=(cov>100);

indels=(cnts(end-1,:,sample)>(cov/10));



for i=1:length(x)
    h=plot(x(i),y(i),'o','MarkerSize', 5, 'MarkerFaceColor', clrsig(psig(i) +1), 'MarkerEdgeColor', clrcov(covsig(i) +1)) ;
    set(h,'ButtonDownFcn',@clicked);
    set(h,'userdata', i) ;
    xlabel(xname)
    ylabel(yname)
    title(SampleNames(sample).Sample)
end

    function clicked(src,event)
        ind = get(src,'userdata') ;
        disp(annotations(ind))
        figure(10);clf;

        %Plot number of calls
        subplot(4,1,1); hold on;
        title(['Number of calls.... ' SampleNames(sample).Sample ' p(no strand bias) = ' num2str(cnts(end-5,ind,sample))]);
        a=squeeze(cnts(1:8,ind,:));
        bar(reshape([a; nan(4,size(a,2))],4,[])','stacked')
        legend('A','T','C','G')
        ylabel('Number of reads')
        
        
        %Plot call quality
        subplot(4,1,2); hold on;
        title(['Average call quality.... ' SampleNames(sample).Sample ' p = ' num2str(cnts(end-4,ind,sample))]);
        a=squeeze(cnts(9:16,ind,:));
        bar(reshape([a; nan(4,size(a,2))],4,[])','grouped')
        legend('Aq','Tq','Cq','Gq')
        ylabel('Average Base Quality')
        
        %Plot mapping quality
        subplot(4,1,3); hold on;
        title(['Average mapping quality.... ' SampleNames(sample).Sample ' p = ' num2str(cnts(end-3,ind,sample))]);
        a=squeeze(cnts(17:24,ind,:));
        bar(reshape([a; nan(4,size(a,2))],4,[])','grouped')
        legend('Am','Tm','Cm','Gm')
        ylabel('Average Mapping Quality')
        
        
        %Plot tail distance
        subplot(4,1,4); hold on;
        title(['Average tail distance.... ' SampleNames(sample).Sample ' p = ' num2str(cnts(end-2,ind,sample))]);
        a=squeeze(cnts(25:32,ind,:));
        bar(reshape([a; nan(4,size(a,2))],4,[])','grouped')
        legend('Atd','Ttd','Ctd','Gtd')
        ylabel('Average Tail Distance')


        
        %Plot MAF in region neighboring locus

        figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/6]);clf;hold on;
       
       % disp(annotations(ind))
        if length(annotations(ind).gene)>0
            protein='';
            for i=1:size(annotations(ind).protein,1)
                protein=[protein ' ' annotations(ind).protein(i,:)];
            end
            title([annotations(ind).gene ' : ' protein '....' num2str(annotations(ind).loc1) '-' num2str(annotations(ind).loc2) ' Strand' num2str(annotations(103).strand+1)])
            plot([annotations(ind).loc1,annotations(ind).loc1],[0,1],strand(2-annotations(ind).strand))
            plot([annotations(ind).loc2,annotations(ind).loc2],[0,1],strand(1+annotations(ind).strand))
                        left=annotations(ind).loc1;
            right=annotations(ind).loc2;
        end
        l=(annotations(ind).pos-window_size):(annotations(ind).pos+window_size);
        for s=2:size(dwindows,3)
            if s~=ind
                plot(l, dwindows(:,ind,s),'.k','MarkerSize', 10) %plot other samples blue
            end
        end
        plot(l,dwindows(:,ind,1),'.g', 'MarkerSize', 15) %plot control in green
        plot(l,dwindows(:,ind,sample),'.m','MarkerSize', 15) %plot sample in magenta
        if (numel(iwindows{ind})) >1
            disp(iwindows{ind}(1,:))
            disp(iwindows{ind}(2,:))
            plot(iwindows{ind}(1,:),iwindows{ind}(2,:),'vk','MarkerFaceColor', 'c','MarkerSize', 25)
        end
        axis([l(1) l(end) .5 1])
        ylabel('MAF')
        
        %Open alignment
        % need to install IGV viewer and set it up for communication in Advanced
        % settings. then remover the comments %*
% 
%         p=pos(:,ind);
%         chr=p(1);
%         position=p(2);
%         
%         t = tcpip('localhost', 60152) ;
%         fopen(t) ;
%         
%         bai_name = [SampleNames(sample).ExperimentFolder '/' SampleNames(sample).Sample '/' SampleNames(sample).AlignmentFolder '/aligned.sorted.bam.bai' ];
%         
%         
%         if ~exist(['../' bai_name ], 'file')
%             error('Create bai files for viewing alignment')
%             eval(['!/opt/bin/samtools index ' bai_name(1:end-4)])
%         end
%         
%         if ~IsGenomeLoaded
%             run_cmd(['genome  Reference_Genomes/' RefGenome '/' RefGenome '.genome' ])
%             IsGenomeLoaded = true ;
%         end
%         
%         run_cmd(['load ' bai_name(1:end-4)]) ;
%         run_cmd(['goto ' ScafNames{chr} ':' num2str(position)])
%         
%     
%         fclose(t);
%         delete(t);
%         clear t
%         
%         function run_cmd(c)
%             disp(c)
%             fprintf(t, c)
%             response = fgetl(t)
%         end
%         
%      
        
    end

end
