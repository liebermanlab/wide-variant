function plot_all_isolate_pairwise_distances(mut_freq, SampleNames)
    % NEEDS TO BE CLEANED UP 
    
    strbit = '-'; 
    all_sites = {}; 
    all_site_names = get_unique_site_names(SampleMames); 
    for n = 1:length(all_site_names)
        all_sites{end+1} = (24*(n-1)+1):24*n; 
    end

    clonal_counts = zeros(length(all_sites));
    dist_matrix = zeros(length(all_sites)); 

    figure; hold on; 
    histbins = [0 1 5:5:50]; 
    % get all within site pairs 
    % cc = hsv(length(all_sites)); 
    cc = distinguishable_colors(length(all_sites)); 

    % perhaps use squareform(pdist(mut_calls'))

    for i = 1:length(all_sites)
        [dist, distm] = calculate_dist_within_site(mut_freq, all_sites{i});
        dist_within_site{i} = dist;
        [n,bins] = hist(dist,histbins); 
        plot(bins,n, 'Color', cc(i,:), 'LineWidth', 2)
        dist_matrix(i,i) = mean(dist); 
        fprintf('Mean distance within site %0.2f\n', mean(dist));     
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
        [n,bins] = hist(dist,histbins); 
    %     plot(bins,n, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '-.')
        dist_matrix(sp(1), sp(2)) = mean(dist);
        fprintf('Mean distance between site %0.2f\n', mean(dist)); 

    end

    legend(all_site_names, 'FontSize', 14); 
    set(gca, 'XTick', 1:10:50, ...
            'FontSize', 20); 


    % ___ Pairwise distances between sites ___ %
    figure; 
    imagesc(dist_matrix);            %# Create a colored plot of the matrix values
    colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                             %#   black and lower values are white)

    textStrings = num2str(dist_matrix(:),'%0.1f');  %# Create strings from the matrix values
    textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
    [x,y] = meshgrid(1:length(all_site_names));   %# Create x and y coordinates for the strings
    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center', ...
                    'FontSize', 16);
    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
    textColors = repmat(dist_matrix(:) > midValue,1,3);  %# Choose white or black for the
                                                 %#   text color of the strings so
                                                 %#   they can be easily seen over
                                                 %#   the background color
    set(hStrings,{'Color'},num2cell(textColors,2));
    set(gca, 'XTick', 1:length(all_site_names), ...
            'YTick', 1:length(all_site_names), ...
            'XTickLabel', all_site_names, ...
            'YTickLabel', all_site_names, ...
            'FontSize', 20)

    
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
end