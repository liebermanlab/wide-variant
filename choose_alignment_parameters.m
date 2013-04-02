function choose_alignment_parameters(filename)


    load(filename);
    
    alignmenttypes=unique({alignments.alignment});
    samples=unique({alignments.sample});
    
    
    percent_aligned=zeros(numel(alignmenttypes),numel(samples));
    
    for i=1:numel(alignments)
        sample=find(strcmp(samples,alignments(i).sample));
        alignment=find(strcmp(alignmenttypes,alignments(i).alignment));
        percent_aligned(alignment,sample)= alignments(i).uniquepairs/alignments(i).totalpairs;
    end
    
    HeatMap(percent_aligned, 'ColumnLabels', samples, 'RowLabels',  alignmenttypes, 'ColumnLabelsRotate', 65)
    set(gca, 'XTick', 1:numel(samples))
    set(gca, 'XTickLabel', samples)
    set(gca, 'YTickLabel', alignmenttypes)
    ylabel(alignmenttypes)

end