function [dist_within_site, dist_between_sites, dist_matrix] = plot_all_isolate_pairwise_distances(mut_freq, SampleNames)
    % USES HAMMING DISTANCE BETWEEN PAIRS, PLOTS MEAN DISTANCE 
  
    unifrac_threshold = 27; % arbitrary 
    
    strbit = '-'; 
    all_sites = {}; 
    all_site_names = get_unique_site_names(SampleNames); 
    for n = 1:length(all_site_names)
        all_sites{end+1} = (24*(n-1)+1):24*n; 
    end

    clonal_counts = zeros(length(all_sites));
    dist_matrix = zeros(length(all_sites)); 
    dist_within_site = cell(length(all_sites),1); 
    
    % Create figure with N subplots for N samples 
%     figure(101); hold on; 
    cc = distinguishable_colors(length(all_sites)); 
    
    
    for i = 1:length(all_sites)
        fh = figure; 
        hold on; 
        
        for sp_i = 1:length(all_sites) 
            
            % ___ PAIRWISE DIST WITH ALL OTHER SITES ___ %
            if sp_i ~= i
                site_i = all_sites{i};
                other_i = all_sites{sp_i};
                dist_other = calculate_dist_between_sites(mut_freq, site_i, other_i); 
            
                % CDF 
                [ycdf, xcdf] = cdfcalc2(dist_other); 
                plot(xcdf, ycdf, 'Color', cc(sp_i,:), 'LineWidth', 1.5, 'LineStyle', ':'); 
            
            % ___ WITHIN SITE ___ %
            else
                [dist_within, distm] = calculate_dist_within_site(mut_freq, all_sites{i});
                dist_within_site{i} = dist_within;
                
                % ADD TO DIST MATRIX
                dist_matrix(i,i) = median(dist_within); 
                
                % CDF
                [ycdf, xcdf] = cdfcalc2(dist_within); 
                plot(xcdf, ycdf, 'Color', cc(i,:), 'LineWidth', 4); 
        
            end
        end
        
        % CONFIGURE SUBPLOT
        hold off; 
        title(all_site_names{i}, 'FontWeight', 'bold', 'FontSize', 16); 
        xlabel('SNP distance between isolates', 'FontWeight', 'bold', 'FontSize', 16); 
        ylabel('Fraction of isolate pairs', 'FontWeight', 'bold', 'FontSize', 16); 
        legend(all_site_names); 
        xlim([0 100]); 
        set(gca, 'XTick', 0:10:100, 'YTick', 0:0.2:1, 'TickLength', [0.005 0.005]); 
        saveas(fh, ['20131025_CDF_' all_site_names{i}], 'png'); 
    end

    % get all between site pairs
    site_pairs = combnk(1:length(all_sites),2); 
    dist_between_sites = cell(length(site_pairs),1); 

    for j = 1:size(site_pairs,1)
        sp = site_pairs(j,:); 
        site1 = all_sites{sp(1)};
        site2 = all_sites{sp(2)}; 
        dist = calculate_dist_between_sites(mut_freq, site1, site2);
        
        dist_between_sites{j} = dist; 
        
        % CDF
%         [ycdf, xcdf] = cdfcalc2(dist); 
%         plot(xcdf, ycdf, 'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-.'); 
        
        % ADD TO DIST MATRIX
        dist_matrix(sp(1), sp(2)) = median(dist);
%         fprintf('Mean distance between site %0.2f\n', mean(dist)); 

    end
    
    % calculate hacked version of unifrac
%     get_hack_unifrac(dist_within_site);
%     get_hack_unifrac(dist_between_sites);
%     
%     % histogram plot details
%     legend(all_site_names, 'FontSize', 14); 
%     xlabel('SNP Distance', 'FontSize', 16, 'FontWeight', 'bold'); 
%     ylabel('# Isolates', 'FontSize', 16, 'FontWeight', 'bold'); 
%     set(gca, 'XTick', 1:10:50, ...
%             'FontSize', 20); 

    % ___ Pairwise distances within/between sites ___ %
%     figure; 
%     imagesc(dist_matrix);            %# Create a colored plot of the matrix values
%     title('Pairwise Distance of Isolates', 'FontWeight', 'bold', 'FontSize', 24); 
%     colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
%                              %#   black and lower values are white)
% 
%     textStrings = num2str(dist_matrix(:),'%0.1f');  %# Create strings from the matrix values
%     textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
%     [x,y] = meshgrid(1:length(all_site_names));   %# Create x and y coordinates for the strings
%     hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
%                     'HorizontalAlignment','center', ...
%                     'FontSize', 16);
%     midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
%     textColors = repmat(dist_matrix(:) > midValue,1,3);  %# Choose white or black for the
%                                                  %#   text color of the strings so
%                                                  %#   they can be easily seen over
%                                                  %#   the background color
%     set(hStrings,{'Color'},num2cell(textColors,2));
%     set(gca, 'XTick', 1:length(all_site_names), ...
%             'YTick', 1:length(all_site_names), ...
%             'XTickLabel', all_site_names, ...
%             'YTickLabel', all_site_names, ...
%             'FontSize', 20)
    
    function uniquekeys = get_unique_site_names(snames)
        allkeys = {}; 
        for u = 1:length(snames)
            sn = snames{u}; 
            keyloc = strfind(sn,strbit);
            skey = sn(1:keyloc-1); 
            allkeys{end+1} = skey;
        end
        uniquekeys = unique(allkeys, 'stable'); % prevents sorting 
    end

    function unifrac_ratios = get_hack_unifrac(dist_cell)
        unifrac_ratios = zeros(length(dist_cell),1); 
        for u = 1:length(dist_cell)
            r = sum(dist_cell{u}>unifrac_threshold)./sum(dist_cell{u}<unifrac_threshold); 
            unifrac_ratio(u) = r; 
        end
    end


end