function observed_muts = div_count_observed_mut_types(annot_mutgenes, type)
    % check that type is either N or S
    if strcmp(type, 'N')
        disp('Getting observed non-synonymous mutations'); 
    elseif strcmp(type, 'S')
        disp('Getting observed synonymous mutations'); 
    else
        disp('Please specify either N or S'); 
        return;   
    end

    observed_muts = zeros(4,4); 
    
    nts = {'A','T','C','G'}; 
    pos = [1, 2, 3, 4]; 
    nt_pos = containers.Map(nts,pos); 
    
    num_muts = size(annot_mutgenes,2); 
    for i = 1:num_muts
        mut_type = annot_mutgenes(i).type; 
        if strcmp(mut_type, 'N')
            % get mutated nucleotide 
            mut_nts = annot_mutgenes(i).nts;
            ref = mut_nts(1);   
            mut = mut_nts(2);  

            % add observed count
            observed_muts(nt_pos(ref), nt_pos(mut)) = observed_muts(nt_pos(ref), nt_pos(mut)) + 1; 
        end
    end
    
end