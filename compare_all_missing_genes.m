function compare_all_missing_genes(gene_coverage, gene_error, genes_nr, ...
                                    allcoverage, SampleNames, allmodes, SortByStd)

    % 2013 by Hattie Chung
    
    scrsz = get(0,'ScreenSize'); 
    
    % make table headers 
    colnames = {'MaxStd', 'Number', 'Gene', 'Start', 'End', 'Annotation', 'Isolate_1', 'Isolate_2'}; 
    nonsamplecols = length(colnames); 
    for s = 1:length(SampleNames)
        colnames{end+1} = SampleNames{s}; 
    end
    
    widths = num2cell([60, 48, 60, 50, 50, 100, 75, 75, ones(1,numel(SampleNames))*60 ]);
    
    % make table
    tabledata = cell(length(genes_nr), length(colnames)); 
    significant_isolates = zeros(length(genes_nr), 3); 
    
    % calculate pairwise STD for all genes
    for genenum = 1:length(genes_nr)
        % Calculate Std 
        std_matrix = compare_single_missing_gene(genenum, gene_coverage, gene_error, genes_nr, ...
                                    allcoverage, SampleNames, allmodes);
        [maxval,ind] = max(std_matrix(:)); 
        [isolate_i, isolate_j] = ind2sub(length(SampleNames), ind); 
        
        % store isolates used for calls
        significant_isolates(genenum,:) = [isolate_i, isolate_j, maxval]; 
        
        % Add All Non-Sample Specific Columns
        tabledata(genenum,1:nonsamplecols) = [  {maxval} ...                        % max std
                                                {genenum} ...
                                                {genes_nr(genenum).gene} ...        % gene name
                                                {genes_nr(genenum).indices(1)} ...  % gene pos
                                                {genes_nr(genenum).indices(2)} ...
                                                {genes_nr(genenum).locustag} ...    % locus 
                                                {SampleNames{isolate_i}} ...        % Isolate 1 used to call
                                                {SampleNames{isolate_j}} ];         % Isolate 2 used to call
        
        % Add coverage at each isolate 
        tabledata( genenum, nonsamplecols+1:nonsamplecols+length(SampleNames) ) = num2cell(gene_coverage( genenum,: )); 
    end
                                
    % sort 
    if SortByStd == 1
        tabledata(isnan(significant_isolates(:,3)),1) = {0}; 
        [sorted_tabledata, sortedpos] = sortrows(tabledata,1); 
        sorted_tabledata = flipdim(sorted_tabledata,1); 
        sortedpos = flipdim(sortedpos,1); 
        significant_isolates_sorted = significant_isolates(sortedpos,:); 
    else
        sorted_tabledata = tabledata; 
        sortedpos = 1:size(tabledata,1); 
    end
    
    % display table
    figure(1001); clf; hold on;
    set(1001, 'Position', [20 110 1250 1200]);
    uicontrol('Style','text','Position',[400 45 120 20],'String','Vertical Exaggeration')
    t = uitable('Units','normalized','Position',[0 0 1 .97], 'Data', sorted_tabledata, ...
        'ColumnName', colnames, ...
        'RowName',[], ...
        'CellSelectionCallback',@missing_clicked, ...
        'ColumnWidth',widths);
    h.checkbox1 = uicontrol('Units','normalized','Style','checkbox','String','Show Alignment in IGV when clicked (must have IGV viewer open already)','Min',0,'Max',1, 'Value',0, 'Position',[0 .97 1 .03]);
    
    % plot histogram of std
    figure(1002); clf; 
    set(1002,'Position', [scrsz(3)*1/2 scrsz(4)*0.3 scrsz(3)/3 scrsz(4)/4]);
    hist(significant_isolates(:,3),60); 
    xlabel('Gene content difference (unit of std)', 'FontWeight', 'bold', 'FontSize', 16); 
    
    % SUBFUNCTION 
    
    function missing_clicked(src, event)
        dt = get(src,'data'); 
        rc = event.Indices ;
        rc(1)
        genenum_i = sortedpos(rc(1)); 
        
        sig_isolates = significant_isolates(genenum_i,:); 
        isol_i = sig_isolates(1); 
        isol_j = sig_isolates(2); 
        max_std = sig_isolates(3); 
        
        % get current position properties
        genepos = genes_nr(genenum_i).indices; 
        if genepos(1)>genepos(2)
            genepos = flipdim(genepos,2);
        end
        if genepos(1) <= 100
            xpos = (genepos(1)):genepos(2)+100; 
            front_offset = 0;
            end_offset = 100; 
        else
            xpos = (genepos(1)-100):(genepos(2)+100); 
            front_offset = 100; 
            end_offset = 100; 
        end
        allmodes_expanded = repmat(allmodes, [1 length(xpos)]); 
        
        
        % plot on fixed figure handle
        figure(121); clf; hold on;
        set(121,'Position',[scrsz(3)*1/2 scrsz(4)*.6 scrsz(3)/3 scrsz(4)/4]);
        
        plot(allcoverage(isol_i, xpos)'./allmodes(isol_i), 'b', 'LineWidth', 1.5); 
        plot(allcoverage(isol_j, xpos)'./allmodes(isol_j), 'k', 'LineWidth', 1.5); 
        legend(SampleNames{isol_i}, SampleNames{isol_j}, ...
                'FontSize', 16, 'FontWeight', 'bold', 'Location', 'NorthWest'); 
        plot((allcoverage(:, xpos)./allmodes_expanded)', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5); 
        % plot again on top of gray lines 
        plot(allcoverage(isol_i, xpos)'./allmodes(isol_i), 'b', 'LineWidth', 1.5); 
        plot(allcoverage(isol_j, xpos)'./allmodes(isol_j), 'k', 'LineWidth', 1.5); 
        
        % plot params
        genelen = genepos(2)-genepos(1); 
        
        % plot gene limits 
        plot([front_offset front_offset], [0 1000], 'r'); 
        plot([genelen+front_offset genelen+front_offset], [0 1000], 'g'); 
        set(gca, 'XTick', 1:500:genelen+front_offset+end_offset, ...
                'XTickLabel', genepos(1)-front_offset:500:genepos(2)+end_offset); 
        xlim([0 genelen+front_offset+end_offset]); 
        ylim([0 5]);
        title(sprintf('Gene %i, Std %0.2f, %s', genenum_i, max_std, genes_nr(genenum_i).product), ...
            'FontWeight', 'bold', 'FontSize', 16); 
        xlabel('Genome Position', 'FontSize', 16, 'FontWeight', 'bold'); 
        ylabel('# Reads/Mode', 'FontSize', 16, 'FontWeight', 'bold'); 

    end

end