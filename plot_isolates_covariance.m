function co_occurrence = plot_isolates_covariance(fixedmut_freq, annotation_genes)
    [num_fixedmuts, num_isolates] = size(fixedmut_freq); 
    co_occurrence = zeros(num_fixedmuts, num_fixedmuts); 
    for iso = 2:num_isolates
        all_pos = find(fixedmut_freq(:,iso)); 
        all_pairwise_pos = combnk(all_pos,2); 
        
        % for all pairwise combinations
        for p = 1:nchoosek(all_pos,2)
            pair = all_pairwise_pos(p,:); 
            co_occurrence(pair(2), pair(1)) = co_occurrence(pair(2),pair(1))+1; 
        end

    end
    
    % get gene names
    mut_genes = annotation_genes.gene; 

    % plot 
    % do some statistics!
    figure; 
    h = imagesc(co_occurrence);
    set(h, 'ButtonDownFcn', @show_selection_info); 
    colormap('summer'); colorbar; 
    set(gca, ...
    'XTick', 1:num_fixedmuts, ... 
    'YTick', 1:num_fixedmuts, ...
    'YTickLabel', mut_genes, ...
    'FontSize', 12, 'FontName', 'Courier'); 

    function show_selection_info(src, event)
        % convert mouse click position to [gene1, gene2]
        point = get(gca, 'CurrentPoint');  
        gene1_pos = floor(point(1,1)+0.5); 
        gene2_pos = floor(point(1,2)+0.5);  
        
        % disp annotation info 
        gene1_pos
        gene2_pos 
        disp(annotation_genes(gene1_pos)); 
        disp(annotation_genes(gene2_pos)); 
    end

end