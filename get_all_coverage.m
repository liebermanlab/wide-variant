function [linear_avg_coverage, log_avg_coverage, zero_cov_positions] = get_all_coverage(dirname, num_expected) 
    % USAGE: get_all_coverage('SA*')
    
    % params 
    threshold = 1; 
    genome_size = 4851126; 
    
    all_dirs = dir(dirname);
    if nargin == 2
        assert(length(all_dirs)==num_expected, 'Incorrect number of samples retrieved'); 
    end
    
    linear_avg_coverage = zeros(length(all_dirs),1);
    log_avg_coverage = zeros(length(all_dirs),1); 
    zero_cov_positions = zeros(length(all_dirs),1);
    
    
    figure(300);
    
    for i = 1:length(all_dirs)    
        sname = all_dirs(i).name; 
        
        % get coverage stats 
        [log_coverage, zero_coverage, log_non_zero, avg_coverage] = get_coverage(sname, threshold); 
        linear_avg_coverage(i) = avg_coverage; 
        zero_cov_positions(i) = length(zero_coverage); 
        log_avg_coverage(i) = mean(log_non_zero); 
        
        % plot coverage 
        figure(300); 
        subplot(6,8,i);
        hist(log_coverage,30); hold on; 
        text(0.5, 8*10^5, sprintf('%0.2f', avg_coverage), ...
            'FontSize', 14, 'FontWeight', 'bold'); 
        text(0.5, 7*10^5, sprintf('%0.2f = %0.2f', mean(log_non_zero), 10^(mean(log_non_zero)) ),...
            'FontSize', 14, 'FontWeight', 'bold'); hold off; 
        title(sname, 'FontSize', 14, 'FontWeight', 'bold'); 
        xlim([0 2.5]); 
        ylim([0 10^6]); 
        
        % print missing coverage
        fprintf('\n\tMissing %0.2f of reference\n', length(zero_coverage)./genome_size); 
        
    end
end