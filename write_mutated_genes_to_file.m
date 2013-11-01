function write_mutated_genes_to_file(annots, outfile)
    % 2013 HC
    
    fid = fopen(outfile, 'w');
    for i = 1:length(annots)
        if strcmp(annots(i).type, 'S') || strcmp(annots(i).type, 'N')
            % only write to file if gene
            fprintf(fid, '%s\t%s\t%s\n', annots(i).type, annots(i).gene, annots(i).protein); 
        end
    end
    fclose(fid); 
end