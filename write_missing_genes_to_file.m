function write_missing_genes_to_file(genes_nr, all_missing, some_missing)
    
    % write uniformly missing genes to file 
    fid1 = fopen('missing_genes_all.csv', 'w');
    fprintf(fid1, '%s,%s,%s,%s\n', 'Start', 'End', 'Gene', 'Locustag'); 
    all_missing_struct = genes_nr(all_missing); 
    for i = 1:numel(all_missing_struct)
        if isempty(all_missing_struct(i).gene)
            gene = '';
        else
            gene = all_missing_struct(i).gene; 
        end
        fprintf(fid1, '%i,%i,%s,%s\n', all_missing_struct(i).loc1, all_missing_struct(i).loc2, ...
                                    gene, all_missing_struct(i).locustag); 
    end
    fclose(fid1);
        
    % write differentially missing genes
    fid2 = fopen('missing_genes_differential.csv', 'w'); 
    fprintf(fid2, '%s,%s,%s,%s\n', 'Start', 'End', 'Gene', 'Locustag'); 
    some_missing_struct = genes_nr(some_missing); 
    for j = 1:numel(some_missing_struct)
        
        if isempty(some_missing_struct(j).gene)
            gene = '';
        else
            gene = some_missing_struct(j).gene; 
        end
        fprintf(fid1, '%i,%i,%s,%s\n', some_missing_struct(j).loc1, some_missing_struct(j).loc2, ...
                                    gene, some_missing_struct(j).locustag); 
    end
    fclose(fid2); 
    
end