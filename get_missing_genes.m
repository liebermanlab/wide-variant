[gene_coverage, gene_names, average_coverages, genes_nr] = calculate_gene_coverage(SampleNames, cds, GenomeLength); 

% only look at where coverage is greater
coverage_threshold = 14; 
gene_coverage_threshold = gene_coverage(:, average_coverages>coverage_threshold); 
samplenames_threshold = SampleNames(average_coverages>coverage_threshold); 

% get genes missing in all strains
genecov_all = sum(gene_coverage_threshold,2);
all_missing = genecov_all==0; 
all_missing_genes = gene_names(all_missing); 

% get genes differentially missing in strains
some_missing = (genecov_all>0) & (genecov_all<length(samplenames_threshold)); 
some_missing_genes = gene_names(some_missing); 

missing_isolates = gene_coverage_threshold(some_missing,:); 

write_missing_genes_to_file(genes_nr, all_missing, some_missing); 
