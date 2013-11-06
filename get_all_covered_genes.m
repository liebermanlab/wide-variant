function [all_coverage_per_bp, cov_modes] = get_all_covered_genes(samplenames, genome_size) 
    % 2013 by Hattie Chung 

    % USAGE: get_all_covered_genes('SA*')
    
    % RETURNS: 
    % ___ all_coverage: [n x GenomeLength] coverage at each bp
    
    % params 
    threshold = 1; 
    all_coverage_per_bp = uint16(zeros(length(samplenames), genome_size));
    cov_modes = zeros(length(samplenames),1); 
    % logical strucuture
    % 1 if position is found, 0 if missing
    for i = 1:length(samplenames)    
        sname = samplenames{i};
        % get covered positions
        fprintf('\nGetting coverage for isolate %i of %i\n', i, length(samplenames))
        [coverage, coverage_mode] = get_covered_genes_sample(sname, threshold); 
        all_coverage_per_bp(i,:) = coverage;
        cov_modes(i) = coverage_mode; 
    end
end