function [linear_avg_coverage, log_avg_coverage, zero_cov_positions] = get_all_coverage(dirname, num_expected) 
    % USAGE: get_all_coverage('SA*')
    
    % params 
    threshold = 1; 
    genome_size = 4851126; 
    
    all_dirs = dir(dirname);
    if nargin == 2
        assert(length(all_dirs)==num_expected, 'Incorrect number of samples retrieved'); 
    end
    all_dirs_sorted = sort_nat({all_dirs.name}); 
    
    linear_avg_coverage = zeros(length(all_dirs),1);
    log_avg_coverage = zeros(length(all_dirs),1); 
    zero_cov_positions = zeros(length(all_dirs),1);
    
    figure(300);
    
    for i = 1:length(all_dirs_sorted)    
        sname = all_dirs_sorted{i}; 
        
        % get coverage stats 
        [log_coverage, zero_coverage, log_non_zero, avg_coverage] = get_coverage(sname, threshold); 
        linear_avg_coverage(i) = avg_coverage; 
        zero_cov_positions(i) = length(zero_coverage); 
        log_avg_coverage(i) = mean(log_non_zero); 
        
        % plot coverage 
        figure(300); 
        subplot(4,6,i);
        hist(log_coverage,30); hold on; 
        text(0.3, 8*10^5, sprintf('%0.2f', avg_coverage), ...
            'FontSize', 10, 'FontWeight', 'bold'); 
        text(0.3, 7*10^5, sprintf('%0.2f = %0.2f', mean(log_non_zero), 10^(mean(log_non_zero)) ),...
            'FontSize', 10, 'FontWeight', 'bold'); hold off; 
        title(sname, 'FontSize', 12, 'FontWeight', 'bold'); 
        xlim([0 2.5]); 
        ylim([0 10^6]); 
        set(gca, 'FontSize', 7, ...
                'XTick', 0:2); 
        set(gcf, 'PaperSize', [12 8], 'PaperPosition', [0 0 12 8]);  
        
        % print missing coverage
        fprintf('\n\tMissing %0.2f of reference\n', length(zero_coverage)./genome_size); 
        
    end
end