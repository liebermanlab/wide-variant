function plot_heatmap_diverse_frequency(mut_freq, mut_pos, annotation_genes, samplenames)
    % sample mutations either fixedmutation or diversemutation 

    Nsample = size(mut_freq,2); 
    
    % get gene names
    mut_genes = {annotation_genes.gene}; 
    no_gene_name = cellfun(@isempty, mut_genes); 
    mut_genes(no_gene_name) = {annotation_genes(no_gene_name).protein}; 
    
    % plot
    figure; 
    imagesc(mut_freq); 
    colormap('summer'); colorbar; 
    set(gca, ...
    'XTick', 1:Nsample, ... 
    'YTick', 1:size(mut_pos,1), ...
    'XTickLabel', samplenames, ... 
    'YTickLabel', mut_genes, ...
    'FontSize', 12, 'FontName', 'Courier'); 
    
end