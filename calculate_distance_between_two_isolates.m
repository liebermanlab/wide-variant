function dist = calculate_distance_between_two_isolates(fixedmutation, isolate1_sample_number, isolate2_sample_number)
    fixedmuts_union = (fixedmutation(:,isolate1_sample_number) | fixedmutation(:,isolate2_sample_number));
    fixedmuts_intersect = (fixedmutation(:,isolate1_sample_number) & fixedmutation(:,isolate2_sample_number)); 
    dist = sum(fixedmuts_union - fixedmuts_intersect); 
end