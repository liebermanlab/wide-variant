% get gene coverage for all samples 
% if exist('gene_coverage_data.mat')
%     fprintf('\nLoading existing gene coverage data...\n'); 
%     load('gene_coverage_data');
% else
%     [gene_coverage, gene_coverage_error, coverage_modes, genes_nr] = calculate_gene_coverage(SampleNames, cds, GenomeLength); 
%     save('gene_coverage_data', 'gene_coverage', 'gene_coverage_error', 'coverage_modes', 'genes_nr'); 
% end


% only take isolate subset that meet coverage threshold
[gene_coverage, gene_coverage_error, coverage_modes, genes_nr, ALLCOVERAGE] = calculate_gene_coverage(SampleNames, cds, GenomeLength); 
coverage_threshold = 10;

gene_coverage_threshold = gene_coverage(:,coverage_modes>coverage_threshold); 
gene_error_threshold = gene_coverage_error(:,coverage_modes>coverage_threshold); 

coverageSampleNames = SampleNames(coverage_modes>coverage_threshold); 
% min_cov_per_gene = min(gene_coverage_threshold,[],2); % need to fix for when min = 0 
% integerize_gene_coverage = gene_coverage_threshold./repmat(min_cov_per_gene, [1 length(coverageSampleNames)]);

%% make clickable table for missing genes 
sort_std = 0; 
compare_all_missing_genes(gene_coverage_threshold, gene_error_threshold, ...
                            genes_nr, ALLCOVERAGE, coverageSampleNames, coverage_modes, sort_std); 

%% plot distribution of gene coverage across isolates
figure, imagesc(gene_coverage_threshold, [0 4]);
xlabel('Isolates', 'FontWeight', 'bold', 'FontSize', 24); 
ylabel('ALL GENES', 'FontWeight', 'bold', 'FontSize', 24); 
title('Not normalized'); 
set(gca, 'YTick', [], ...
        'FontSize', 14, ...
        'TickLength', [0.001 0.001]);

%% normalizing by GENE coverage mode 
m_all = zeros(size(gene_coverage_threshold,2),1); 
for i = 1:size(gene_coverage_threshold,2)
    [n,bin] = hist(gene_coverage_threshold(:,i),[0:0.1:10]); 
    m_i = bin(n==max(n(5:end))); 
    m_all(i) = m_i; 
end

test_norm = gene_coverage_threshold./repmat(m_all,[1 size(gene_coverage_threshold,1)])'; 
figure, imagesc(test_norm, [0 4]); 
% colormap('bone'); 
xlabel('Isolates', 'FontWeight', 'bold', 'FontSize', 24); 
ylabel('ALL GENES', 'FontWeight', 'bold', 'FontSize', 24); 
title('Normalize by GENE COV MODE'); 
set(gca, 'YTick', [], ...
        'FontSize', 14, ...
        'TickLength', [0.001 0.001]);    
    
