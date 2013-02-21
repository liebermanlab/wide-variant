
function div_clickable_scatter(x,y, xname, yname, sample, toplot, params, covthresholds, cnts, dwindows, iwindows, cwindows, pos, annotations, RefGenome, ScafNames, ChrStarts, SampleNames)

scrsz = get(0,'ScreenSize');

IsGenomeLoaded = false ;
window_size=floor(size(dwindows,1)/2);


tamismagenta=[.85 .31 .85];
tamisgreen=[145 220 83]/256;


%Calculate important things, used later for plotting window
[dmaf, ~, ~] = div_major_allele_freq(cnts);
d=div_test_thresholds(cnts,params, covthresholds);
goodmaf=zeros(size(d)); goodmaf(d>0)=dmaf(d>0);
p=chrpos2index(pos',ChrStarts);

good=d(:,sample)>0;


fign=1000+floor(100*rand());
figure(fign); clf ; hold on;
set(fign,'Position',[scrsz(3)/20 scrsz(4)/8 scrsz(3)/5 scrsz(3)/5]);clf;hold on;



fixed=zeros(size(toplot));

%add randomization 
for i=1:numel(toplot)
    if x(i)> .9 &  y(i)>.9
        x(i)=x(i)-(rand()*.03)+.01;
        y(i)=y(i)-(rand()*.04)+.02;
        fixed(i)=1;
    end
end

x=x-(rand(size(x))*.01)+.01;

p1=plot(x(toplot & good & ~fixed),y(toplot & good & ~fixed),'o');
set(p1,'MarkerSize', 4, 'MarkerFaceColor', tamismagenta, 'MarkerEdgeColor', tamismagenta, 'ButtonDownFcn',@clicked) ;
p2=plot(x(~good & toplot& ~fixed),y(~good & toplot & ~fixed),'s');
set(p2,'MarkerSize', 4, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', rgb('Gray'),'ButtonDownFcn',@clicked);

p3=plot(x(~good & toplot &fixed),y(~good & toplot & fixed),'o');
%set(p3,'MarkerSize', 4, 'MarkerFaceColor', rgb('Magenta'), 'MarkerEdgeColor',rgb('Magenta')),'ButtonDownFcn',@clicked);
set(p3,'MarkerSize', 4, 'MarkerFaceColor', tamisgreen, 'MarkerEdgeColor',tamisgreen,'ButtonDownFcn',@clicked);

disp(sum(~good & toplot))



search_ind = find(toplot) ;

plot([-.5 1.5], [params.minorfreqthreshold params.minorfreqthreshold], '--', 'Color', rgb('Black'))
%plot([-.5 1.5], [.025 .025], '--', 'Color', rgb('Black'))

axis([-.02 1.02 -0.01 1.03])

xlabel(xname)
ylabel(yname)
set(gca,'Xtick',0:.1:1)
set(gca,'Ytick',0:.1:1)
grid on
set(gca, 'GridLineStyle', ':');


title(SampleNames(sample).Sample)



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
        [~,ind] = min((x(search_ind)-x0).^2+(y(search_ind)-y0).^2) ;
        
        ind = search_ind(ind) ;
        
        disp(ind);
        
        disp(annotations(ind))
        
        chr=annotations(ind).scafold;
        position= annotations(ind).pos;
        div_display_thresholds(cnts(:,ind,sample), params, dmaf(ind,1),covthresholds)
        disp(cnts(:,ind,sample))
        
        %Bar charts of counts
        div_bar_charts(squeeze(cnts(:,ind,:)), sample, {SampleNames.Sample})
     
        
        
        %Plot MAF in region neighboring locus
        region=(find(p>p(ind)-window_size,1):find(p<p(ind)+window_size,1,'last'));
      %  iwindows{ind}(2)
        div_maf_window(annotations(ind), p(ind)-ChrStarts(chr), window_size, iwindows{ind}, squeeze(dwindows(:,ind,:)), p(region)-ChrStarts(chr), goodmaf(region,:), {SampleNames.Sample})
        div_cov_window(annotations(ind), p(ind)-ChrStarts(chr), window_size, squeeze(cwindows(:,ind,:)), {SampleNames.Sample})

        
        
        
        %Open alignment
        % need to install IGV viewer and set it up for communication in Advanced
        %settings. then remover the comments %*
        
        show_alignment=0;
        if show_alignment
            ptemp=pos(ind,:);
            chr=ptemp(1);
            position=ptemp(2);
            
            t = tcpip('localhost', 60152) ;
            fopen(t) ;
            
            bai_name = [SampleNames(sample).ExperimentFolder '/' SampleNames(sample).Sample '/' SampleNames(sample).AlignmentFolder '/aligned.sorted.bam.bai' ];
            
            
            if ~exist(['../' bai_name ], 'file')
                error('Create bai files for viewing alignment')
                eval(['!/opt/bin/samtools index ' bai_name(1:end-4)])
            end
            
            if ~IsGenomeLoaded
                run_cmd(['genome  Reference_Genomes/' RefGenome '/' RefGenome '.genome' ])
                run_cmd(['load ' bai_name(1:end-4)]) ;
                IsGenomeLoaded = true ;
            end
            
            
            run_cmd(['goto ' ScafNames{chr} ':' num2str(position)])
            
            
            fclose(t);
            delete(t);
            clear t
            
        end
        
        function run_cmd(c)
            disp(c)
            fprintf(t, c)
            response = fgetl(t)
        end
        
        
        
        
    end

end
