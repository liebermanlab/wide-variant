function pairwise_std = compare_single_missing_gene(genenum, gene_coverage, gene_error, genes_nr, ...
                                                        allcoverage, SampleNames, allmodes)

    cov = gene_coverage(genenum,:); 
    err = gene_error(genenum,:); 
    pairwise_std = zeros(length(cov)); 
    
    for i = 1:length(cov)
        for j = 1:length(cov) 
            pairwise_std(i,j) = (cov(i)-cov(j))./sqrt(err(i).^2+err(j).^2); 
        end
    end
    
    [maxval,ind] = max(pairwise_std(:)); 
    [colony_i, colony_j] = ind2sub(length(SampleNames),ind);

    
    
end