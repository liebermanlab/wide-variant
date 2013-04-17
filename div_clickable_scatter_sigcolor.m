
function div_clickable_scatter_sigcolor(x,y, xname, yname, sample, params, covthresholds, cnts, fwindows, cwindows, pos, annotations, RefGenome, ScafNames, ChrStarts, SampleInfo)




IsGenomeLoaded = false ;
window_size=floor(size(fwindows,1)/2);


%Calculate important things, used later for plotting window
d=div_test_thresholds(cnts,params, covthresholds);
good=d(:,sample)>0;
p=chrpos2index(pos',ChrStarts);
[dmaf, ~, ~] = div_major_allele_freq(cnts);
goodmaf=zeros(size(d)); goodmaf(d>0)=dmaf(d>0);

figure;
clf ; 
subplot(2,1,1); hold on; 

search_ind = 1:numel(good) ;
%search_ind = find(good) ;
%search_ind = find(~good) ;

plot(x(~good),y(~good),'o','MarkerSize', 4, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'ButtonDownFcn',@clicked) ;

plot(x(good),y(good),'o','MarkerSize', 4, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'ButtonDownFcn',@clicked) ;




axis([.5 1 .5 1])

xlabel(xname)
ylabel(yname)
title(SampleInfo(sample).Sample)


% Decide where to click

% Create the button group.
h = uibuttongroup('visible','off','Position',[0 0 .5 .4]);
% Create three radio buttons in the button group.
u0 = uicontrol('Style','Radio','String','Inspect nearest position regardless of color',...
    'pos',[10 75 1100 30],'parent',h,'HandleVisibility','off');
u1 = uicontrol('Style','Radio','String','Inspect nearest called position (red)',...
    'pos',[10 50 1100 30],'parent',h,'HandleVisibility','off');
u2 = uicontrol('Style','Radio','String','Inspect nearest uncalled position (black)',...
    'pos',[10 25 1100 30],'parent',h,'HandleVisibility','off');
% Initialize some button group properties. 
set(h,'SelectionChangeFcn',@selcbk);
set(h,'SelectedObject',[]);  % No selection
set(h,'Visible','on');


    function selcbk(source,eventdata)
        
        if strcmp(get(get(source,'SelectedObject'),'String'),'Inspect nearest called position (red)');
            search_ind = find(good);
        elseif strcmp(get(get(source,'SelectedObject'),'String'),'Inspect nearest uncalled position (black)');
            search_ind = find(~good);
        end
    end

    function clicked(src,event)
        
        
        search_ind = find(~good) ;

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
        disp(annotations(ind))
        div_display_thresholds(cnts(:,ind,sample), params, dmaf(ind,1), covthresholds)
        disp(cnts(:,ind,sample))


        %Bar charts of counts
        div_bar_charts(squeeze(cnts(:,ind,:)), sample, {SampleInfo.Sample})
     
        
        
        %Plot MAF in region neighboring locus
        region=(find(p>p(ind)-window_size,1):find(p<p(ind)+window_size,1,'last'));
        if numel(fwindows)>1
           
            div_maf_window(annotations(ind), p(ind)-ChrStarts(chr), window_size, [], ...
               squeeze(fwindows(:,ind,:)), p(region)-ChrStarts(chr), goodmaf(region,:), ...
               {SampleInfo.Sample},sample,1) % showlegends)
           
           div_cov_window(annotations(ind), p(ind)-ChrStarts(chr), window_size, ...
               squeeze(cwindows(:,ind,:)), {SampleInfo.Sample},sample, 1) %showlegends)
        end
        

        
        
        %Open alignment
        % need to install IGV viewer and set it up for communication in Advanced
        %settings. then remover the comments %*
        
        
        if show_alignment==1
            
            t = tcpip('localhost', 60152) ;
            fopen(t) ;
            
            bai_name = [SampleInfo(sample).ExperimentFolder '/' SampleInfo(sample).Sample '/' SampleInfo(sample).AlignmentFolder '/aligned.sorted.bam.bai' ];
            
            
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
