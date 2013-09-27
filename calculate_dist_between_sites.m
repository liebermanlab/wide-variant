function all_distances = calculate_dist_between_sites(fixedmuts, site1_isolates, site2_isolates) 
    % HC 7/26/2013
    % calculates pairwise distance between isolates in two different sites
    
    % calculate all combinations
    both_sites = combnk([site1_isolates site2_isolates],2);
    % calculate all combinations in site 1
    site1 = combnk(site1_isolates,2); 
    % calculate all combinations in site 2
    site2 = combnk(site2_isolates,2); 
    
    % use setdiff to remove site1 + site2 from all_sites
    between_sites = setdiff(setdiff(both_sites, site1, 'rows'), site2, 'rows');
    num_pairs = size(between_sites,1); 
    % calculate all pairwise dist
    all_distances = zeros(num_pairs,1); 
    for i = 1:num_pairs
        pair = between_sites(i,:);
        dist = calculate_distance_between_two_isolates(fixedmuts, pair(1), pair(2)); 
        all_distances(i) = dist; 
    end
end