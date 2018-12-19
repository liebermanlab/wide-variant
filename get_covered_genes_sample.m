function [coverage, coverage_mode] = get_covered_genes_sample(sample_name, threshold, genome_size)
    % USAGE: get_coverage('SA11-1') 

    if nargin < 2
        threshold = 1; 
    end
    
    
    % load diversity.mat for each file 
    clear diversity_file diversity 

    fprintf('\nLoading file for %s\n', sample_name); 
    diversity_file = strcat(sample_name, '/diversity.mat'); 
    if exist(diversity_file)
        diversity = load(diversity_file);
        countsdata = diversity.data;
    
        % calculate coverage
        coverage = sum(countsdata(1:8,:));

        average_coverage = mean(coverage); 

        % get mode
        modebins = [0:1:max(coverage)]; 
        [n,bins] = hist(coverage,modebins);
        coverage_mode = bins(n==max(n(5:end))); 
    else
        coverage = uint16(1,zeros(genome_size)); 
        coverage_mode = 0; 
    end
   
end