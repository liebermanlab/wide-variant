function [all_distances, dist_matrix] = calculate_dist_within_site(fixedmuts, site_isolates)
    % HC 7/26/2013
    % given site_isolates, an array with positions for isolate #s in
    % mut_freq, calculates distribution of pairwise distances
    
    % calculate all pairwise combinations
    all_pairs = combnk(site_isolates,2); 
    num_pairs = size(all_pairs,1); 
    
    % hack to get position for dist_matrix
    all_pairs_index = combnk(1:length(site_isolates),2); 
    
    % get distances
    all_distances = zeros(num_pairs,1);
    dist_matrix = zeros(length(site_isolates),length(site_isolates)); 
    for i = 1:num_pairs
        pair = all_pairs(i,:); 
        pair_index = all_pairs_index(i,:); 
        dist = calculate_distance_between_two_isolates(fixedmuts, pair(1), pair(2)); 
        dist_matrix(pair_index(1), pair_index(2)) = dist; 
        all_distances(i) = dist;
    end
    
end