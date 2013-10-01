function [coverage, average_coverage] = get_covered_genes_sample(sample_name, threshold)
    % USAGE: get_coverage('SA11-1') 

    if nargin < 2
        threshold = 1; 
    end
    
    % load diversity.mat for each file 
    clear diversity_file diversity 

    fprintf('\nLoading file for %s\n', sample_name); 
    diversity_file = strcat(sample_name, '/diversity.mat'); 
    diversity = load(diversity_file);
    counts = diversity.data;
    
    % calculate coverage
    linear_coverage = sum(counts(1:8,:));
    
    % average coverage
    coverage = linear_coverage>threshold;
    
    average_coverage = mean(linear_coverage(linear_coverage>threshold)); 
end