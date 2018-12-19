function [multimut_genenums, genotypes] = get_multiple_mutated_loci(genes_mut_cnt, genes, unique_genes, Callsgood)
    % 2013 Hattie Chung 
    
    multimut_genenums = unique_genes(genes_mut_cnt(:,2)>1); 
    multimut_positions = cell(length(multimut_genenums),1);
    genotypes = cell(length(multimut_genenums),1); 
    genotype_degenerates = cell(length(multimut_genenums),1); 

    for m = 1:length(multimut_genenums);
        cur_gene_dict = containers.Map(); 
        degenerates = []; 

        % get all structure entries that match gene number 
        mutpositions = find(genes==multimut_genenums(m));
        multimut_positions{m} = mutpositions; 

        % ALL "genotypes" 
        mutcalls_m = Callsgood(mutpositions,:);

        for x = 1:size(mutcalls_m,2)
            g_x = mutcalls_m(:,x); 

            % add to key if no N's and not in current hash table 
            if ~isKey(cur_gene_dict,g_x) && isempty(strfind(g_x','N'))
                cur_gene_dict(g_x) = [x]; 
            elseif ~isKey(cur_gene_dict,g_x)
                degenerates(end+1) = [x];
            else
                cur_gene_dict(g_x) = [cur_gene_dict(g_x) x]; 
            end
        end

        genotypes{m} = cur_gene_dict; 
        genotype_degenerates{m} = degenerates; 
        
    end
end