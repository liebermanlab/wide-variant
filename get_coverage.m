function [log_coverage, zero_coverage, non_zero_coverage, average_coverage] = get_coverage(sample_name, threshold)
    % USAGE: get_coverage('SA11-1') 

    if nargin < 2
        threshold = 1; 
    end
    log_threshold = log10(threshold); 
    
    % load diversity.mat for each file 
    clear diversity_file diversity 

    fprintf('\nLoading file for %s\n', sample_name); 
    diversity_file = strcat(sample_name, '/diversity.mat'); 
    diversity = load(diversity_file); 
    counts = diversity.data;
    
    % calculate coverage
    linear_coverage = sum(counts(1:8,:)); 
    log_coverage = log10(linear_coverage+1); 
    
    % average coverage
    zero_coverage = find(log_coverage<=log_threshold); 
%     zero_coverage = find(coverage<=threshold); 
    
    % print average coverage
    non_zero_coverage = log_coverage(log_coverage>log_threshold); 

    average_coverage = mean(linear_coverage(linear_coverage>threshold)); 
    fprintf('\n\tSample %s average %0.2f\n', sample_name, average_coverage); 
    fprintf('\n\tSample %s log average %0.2f\n', sample_name, mean(non_zero_coverage)); 
end