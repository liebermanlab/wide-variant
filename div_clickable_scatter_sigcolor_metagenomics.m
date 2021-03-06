
function div_clickable_scatter_sigcolor_metagenomics(x,y, xname, yname, sample, params, cnts, windows, pos, ChrStarts, SampleNames) %


IsGenomeLoaded = false ;
window_size=floor(size(windows,1)/2);


%Calculate important things, used later for plotting window
d=div_test_thresholds_metagenomics(cnts,params);
good=d(:,sample)>0;
p=chrpos2index(pos',ChrStarts);
[dmaf, ~, ~] = div_major_allele_freq(cnts);
goodmaf=zeros(size(d)); goodmaf(d>0)=dmaf(d>0);

figure;
clf ; hold on

search_ind = find(good) ;

plot(x(~good),y(~good),'o','MarkerSize', 3, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'ButtonDownFcn',@clicked) ;

plot(x(good),y(good),'.','MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'ButtonDownFcn',@clicked) ;



axis([.5 1 .5 1])

xlabel(xname)
ylabel(yname)
title(SampleNames(sample).Sample)

    function clicked(src,event)
        
        
        
        show_alignment=0

        
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
        [~,ind] = min((x(search_ind)-x0).^2+(y(search_ind)-y0).^2) ;
        
        ind = search_ind(ind) ;
        ptemp=pos(ind,:);
        chr=ptemp(1);
        position=ptemp(2);
        
        
        
        disp(ind)
        disp(cnts(:,ind,sample))


        %Bar charts of counts
        div_bar_charts(squeeze(cnts(:,ind,:)), sample, SampleNames(sample).Sample)
     
        
        
        %Plot MAF in region neighboring locus
        
        region=(find(p>ChrStarts(chr),1):find(p<(ChrStarts(chr+1)-1),1,'last'));
        div_maf_window_metagenomics(chr, position, squeeze(windows(:,ind,:)), p(region)-ChrStarts(marker), goodmaf(region,:), {SampleNames.Sample})
        

        
        
        %Open alignment
        % need to install IGV viewer and set it up for communication in Advanced
        %settings. then remover the comments %*
        
        
%         if show_alignment
%             
%             t = tcpip('localhost', 60152) ;
%             fopen(t) ;
%             
%             bai_name = [SampleNames(sample).ExperimentFolder '/' SampleNames(sample).Sample '/' SampleNames(sample).AlignmentFolder '/aligned.sorted.bam.bai' ];
%             
%             
%             if ~exist(['../' bai_name ], 'file')
%                 error('Create bai files for viewing alignment')
%                 eval(['!/opt/bin/samtools index ' bai_name(1:end-4)])
%             end
%             
%             if ~IsGenomeLoaded
%                 run_cmd(['genome  Reference_Genomes/' RefGenome '/' RefGenome '.genome' ])
%                 run_cmd(['load ' bai_name(1:end-4)]) ;
%                 IsGenomeLoaded = true ;
%             end
%             
%             
%             run_cmd(['goto ' ScafNames{chr} ':' num2str(position)])
%             
%             
%             fclose(t);
%             delete(t);
%             clear t
%             
%         end
%         
%         function run_cmd(c)
%             disp(c)
%             fprintf(t, c)
%             response = fgetl(t)
%         end
%         
%         
%         
%         
%     end

end
