% each site's isolates
% SA11 = 1:24; 
% SA12 = 25:48; 
% SA13 = 49:72;
% SA15 = 73:96; 
% SA18 = 97:120;
% SA19 = 121:126; 
% all_sites = {SA11, SA12, SA13, SA15, SA18, SA19}; 
% dist_within_site = cell(length(all_sites),1); 

all_sites = {}; 
all_site_names = {'SA1', 'SA4', 'SA5', 'SA7', 'SA8', 'SA9', 'SA10', 'SA19', 'SA20', 'SA13', ...
            'SA15', 'SA18', 'SA11', 'SA12'}; 
for n = 1:length(all_site_names)
    all_sites{end+1} = (24*(n-1)+1):24*n; 
end


dist_matrix = zeros(length(all_sites)); 
% get all within site pairs 
cc = hsv(length(all_sites)); 
for i = 1:length(all_sites)
    dist = calculate_dist_within_site(mut_freq, all_sites{i});
    dist_within_site{i} = dist; 
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
    dist_matrix(sp(1), sp(2)) = mean(dist); 
    fprintf('Mean distance between site %0.2f\n', mean(dist)); 
    
end

figure; 
imagesc(dist_matrix);            %# Create a colored plot of the matrix values
colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                         %#   black and lower values are white)

textStrings = num2str(dist_matrix(:),'%0.1f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:14);   %# Create x and y coordinates for the strings
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