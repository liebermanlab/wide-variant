function [gene_cov_all_isolates, gene_err_all_isolates, cov_modes, cds_nonredundant, coverage_per_bp] = calculate_gene_coverage(samplenames, cds, genomesize)
    % 2013 by Hattie Chung

    % params
    readlen = 100; 

    % initialize 
    cds_nonredundant_names = cell(0,1); 
    gene_cov_all_isolates = uint16(zeros(0,length(samplenames)));
    gene_err_all_isolates = zeros(0,length(samplenames)); 
    
    % load coverage for isolates
    [coverage_per_bp, cov_modes] = get_all_covered_genes(samplenames, genomesize); 
    
    % annotate coding positions of chromosome
    coding_sequences = cds{1}; 
    cds_nonredundant = coding_sequences; 
    
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
            sum_of_calls_in_gene = coverage_per_bp(:,gene_begin:gene_end); 
            reads_per_gene = (1./readlen)*sum(sum_of_calls_in_gene,2);
            
            % normalize 
            coverage_per_gene = (readlen./gene_len)*(reads_per_gene./cov_modes); 
            gene_cov_all_isolates(end+1,:) = coverage_per_gene; 
            
            error_per_gene = (readlen./gene_len)*(sqrt(reads_per_gene)./cov_modes); 
            gene_err_all_isolates(end+1,:) = error_per_gene; 
            
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