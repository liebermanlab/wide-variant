function all_distances = calculate_dist_within_site(fixedmuts, site_isolates)
    % HC 7/26/2013
    % given site_isolates, an array with positions for isolate #s in
    % mut_freq, calculates distribution of pairwise distances
    
    % calculate all papirwise combinations
    all_pairs = combnk(site_isolates,2); 
    num_pairs = size(all_pairs,1); 
    % get distances
    all_distances = zeros(num_pairs,1); 
    for i = 1:num_pairs
        pair = all_pairs(i,:); 
        dist = calculate_distance_between_two_isolates(fixedmuts, pair(1), pair(2)); 
        all_distances(i) = dist; 
    end
    
end