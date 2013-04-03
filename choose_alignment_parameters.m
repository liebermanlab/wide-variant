function choose_alignment_parameters(directory)


    load(directory);
    
    alignmenttypes=unique({alignments.alignment});
    samples=unique({alignments.sample});
    
    
    percent_aligned=zeros(numel(alignmenttypes),numel(samples));
    
    for i=1:numel(alignments)
        sample=find(strcmp(samples,alignments(i).sample));
        alignment=find(strcmp(alignmenttypes,alignments(i).alignment));
        percent_aligned(alignment,sample)= alignments(i).uniquepairs/alignments(i).totalpairs;
    end
    
    % use imagesc instead of HeatMap
    figure; 
    lims = [0 1]; 
    h = imagesc(percent_aligned, lims);
    colorbar; 
    set(gca, 'XTickLabel', samples); 
    set(gca, 'XTick', 1:numel(samples)); 
    set(gca, 'YTickLabel', []); 
    set(gca, 'YTick', []); 
    th = title(alignmenttypes, 'FontSize', 24, 'FontWeight', 'bold'); 
    set(th, 'Interpreter', 'none'); 
end