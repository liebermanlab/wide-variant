function [cds_isolates, cds_nonredundant_names, avg_cov, cds_nonredundant] = calculate_gene_coverage(samplenames, cds, genomesize)
    cds_nonredundant_names = cell(0,1); 
    cds_isolates = zeros(0,length(samplenames)); 
    arbitrary_threshold = 0.9; 
    
    % load coverage for isolates
    [coverage_positions, avg_cov] = get_all_covered_genes(samplenames, genomesize); 
    
    % annotate coding positions of chromosome
    coding_sequences = cds{1}; 
    cds_nonredundant = coding_sequences; 
    % initialize 

    old_locustag = ''; 
    for c = 1:length(coding_sequences) 
        pos = coding_sequences(c).indices; 
        locustag = coding_sequences(c).locustag; 
        gene_begin = min(pos); 
        gene_end = max(pos); 
        gene_len = gene_end - gene_begin; 
        gene_name = coding_sequences(c).gene;
        
        % add CDS entry, if non-redundant
        if ~strcmp(locustag, old_locustag)
            cur_gene_coverage = sum(coverage_positions(:,gene_begin:gene_end),2);
            isolates_gene_found = cur_gene_coverage > arbitrary_threshold*gene_len;
            cds_isolates(end+1,:) = isolates_gene_found; 
            
            if ~isempty(gene_name)
                cds_nonredundant_names{end+1} = gene_name; 
            else
                cds_nonredundant_names{end+1} = locustag; 
            end
        else
            cds_nonredundant(c) = []; 
        end
        
        old_locustag = locustag; 
    end
    
end