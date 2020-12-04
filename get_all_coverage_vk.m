function [all_coverage_per_bp, cov_modes, all_maf_per_bp] = get_all_coverage_vk(samplefilenames, genome_size) 
    % samplefilenames -- cell array of paths to diversity mat files  
    %OUTPUT: where n is number of samples 
    % all_coverage_per_bp: [n x GenomeLength] coverage at each bp
    % cov_modes: [n x 1] coverage mode for each sample  
    
    % params 
    all_coverage_per_bp = uint16(zeros(length(samplefilenames), genome_size));
    all_maf_per_bp = uint16(zeros(length(samplefilenames), genome_size));
    cov_modes = zeros(length(samplefilenames),1); 
    % logical strucuture
    % 1 if position is found, 0 if missing
    for i = 1:length(samplefilenames)    
        filename = samplefilenames{i};
        % get covered positions
        fprintf('\nGetting coverage for isolate %i of %i\n', i, length(samplefilenames))
        [coverage, coverage_mode, maf] = get_sample_coverage(filename, genome_size);  %fprintf(1,['\n' num2str(numel(coverage)) ' ' num2str(genome_size)]);
        all_coverage_per_bp(i,:) = coverage;
        all_maf_per_bp(i,:)=1000*maf;
        cov_modes(i) = coverage_mode(1); 
    end
end