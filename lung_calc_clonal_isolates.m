function [isolate_type, uniquetypecounts, unique_genotypes]  = lung_calc_clonal_isolates(mut_freq, threshold)
    
    % SNP distance threshold for "clonality" 
    if nargin < 2
        threshold = 3; 
    end

    % Find number of unique genotypes in each site
    siteisolates = mut_freq; 
    [~, num_isolates] = size(siteisolates); 
    isolate_type = zeros(num_isolates,1); 
    for i = 1:num_isolates
        genotype_i = siteisolates(:,i);
        for j = 1:num_isolates
            genotype_j = siteisolates(:,j); 
            dist_i = abs(sum(genotype_i-genotype_j));
            
            % compare two genotypes
            if dist_i < threshold
                if isolate_type(min(i,j)) ~= 0
                    isolate_type([i j]) = isolate_type(min(i,j));
                else
                    isolate_type([i j]) = min(i,j); 
                end
            end
        end
    end
    
    % calculate number of isolate in each unique type 
    uniquetypes = unique(isolate_type);
    uniquetypecounts = zeros(length(uniquetypes),1); 
    for u = 1:length(unique(isolate_type))
        this_type_num = length(find(isolate_type==uniquetypes(u)));
        uniquetypecounts(u) = this_type_num; 
    end
    
    % extract unique genotypes
    unique_genotypes = siteisolates(:,uniquetypes); 
    
end