function [covered_positions, avg_coverages] = get_all_covered_genes(samplenames, genome_size) 
    % USAGE: get_all_coverage('SA*')
    
    % params 
    threshold = 1; 
    covered_positions = repmat(uint8(0), length(samplenames), genome_size);
    avg_coverages = zeros(length(samplenames),1); 
    % logical strucuture
    % 1 if position is found, 0 if missing
    for i = 1:length(samplenames)    
        sname = samplenames{i};
        % get covered positions
        fprintf('\nCalculating missing positions for isolate %i of %i\n', i, length(samplenames))
        [coverage, average_coverage] = get_covered_genes_sample(sname, threshold); 
        covered_positions(i,:) = coverage;
        avg_coverages(i) = average_coverage; 
    end
end