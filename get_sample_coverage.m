function [coverage, coverage_mode, MAF] = get_sample_coverage(filename,genome_size) %sceond input not needed but kept for backwards compatibility

    
    diversity_file = strcat(filename); 
    
    if exist(diversity_file, 'file')
        diversity = load(diversity_file);
        countsdata = diversity.data;
    
        % calculate coverage
        coverage = sum(countsdata(1:8,:));
        average_coverage = mean(coverage); 
        coverage_mode=mode(coverage);
        
        
        % calculate major allele frequency
        c=countsdata(1:4,:)+countsdata(5:8,:);
        [sorted, sortedpositions] = sort(c,1);
        maxcount = sorted(end,:,:);
        MAF=double(maxcount)./sum(c,1);
        MAF(isnan(MAF))=0; %set to 0 to indicate no data

    else
        error(['No diversity.mat file found at :' diversity_file])
    end
   
end