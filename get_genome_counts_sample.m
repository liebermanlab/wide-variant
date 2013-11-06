function [genome_calls] = get_genome_counts_sample(sample_name, threshold)
    % USAGE: get_coverage('SA11-1') 

    if nargin < 2
        threshold = 1; 
    end
    
    
    % load diversity.mat for each file 
    clear diversity_file diversity 

    fprintf('\nLoading file for %s\n', sample_name); 
    diversity_file = strcat(sample_name, '/diversity.mat'); 
    diversity = load(diversity_file);
    countsdata = diversity.data;
    
    % calculate coverage
    genome_calls = countsdata(1:4,:)+countsdata(5:8,:);
%     genome_calls_sum = sum(genome_calls,1); 
%     average_coverage = mean(genome_calls_sum); 
    
    % get mode
%     modebins = [0:1:max(genome_calls_sum)]; 
%     [n,bins] = hist(genome_calls_sum,modebins);
%     coverage_mode = bins(n==max(n(5:end))); 
   
end