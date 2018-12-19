function get_isolates_with_single_loci_mut()
    % get loci numbers that have one mut
    single_mut_loci_num = genes_mutation_count(genes_mutation_count(:,2)==1,1); 
    % compare to all loci numbers
    all_gene_nums = [annotation_all.gene_num];
    [~, ~, i_allgene] = intersect(single_mut_loci_num, all_gene_nums); 
    % extract annotation structure for single mut loci
    annotation_singlemut = annotation_all(i_allgene); 

    % get reference at these positions 
    ancnt_good = ancnt(MutQual>qual_0); 
    ancnt_singlemut = ancnt_good(i_allgene); 
    single_genotypes = cell(length(i_allgene),1); 

    % desize Callsgood for only single mut loci 
    Callsgood_single = Callsgood(i_allgene, :);

    for s = 1:length(annotation_singlemut)
        mutcalls_m = Callsgood(s,:); 

        for x = 1:size(mutcalls_m,2)
            g_x = mutcalls_m(s,x); 
            anc_x = ancnt_singlemut(s); 
            if ~strcmp(g_x, anc_x) 
                % blah 
            end
        end
    end
end