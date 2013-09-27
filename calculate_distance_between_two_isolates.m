function dist = calculate_distance_between_two_isolates(fixedmutation, isolate1, isolate2)
    fixedmuts_union = (fixedmutation(:,isolate1) | fixedmutation(:,isolate2));
    fixedmuts_intersect = (fixedmutation(:,isolate1) & fixedmutation(:,isolate2)); 
    dist = sum(fixedmuts_union - fixedmuts_intersect); 
end