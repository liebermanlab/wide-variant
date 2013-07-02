function [mut_freq, mut_pos, annotation_genes] = get_called_mutations(samplemutations, annotation_all, mutAF, samplenames) 
    
    Nsample = size(samplemutations,2); 
    
    % get diverse positions
    mut_pos = []; 
    for sample = 1:Nsample
        positions = find(samplemutations(:,sample)~=0);
        mut_pos = union(mut_pos, positions); 
    end
    
    % get diverse freq
    mut_freq = mutAF(mut_pos,:); 
    
    % get gene names
    annotation_genes = annotation_all(mut_pos); 
    
end