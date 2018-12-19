% Meant to run from inside the case folder, such that sample_i/diversity.mat exists for all samples
% This is overkill, but need SampleNames and GenomeLength. Can also pull from sample_names.csv
run_postfix='14_06';
%load(['mutation_table_' run_postfix])


SampleInfo = read_sample_names ;
SampleNames={SampleInfo(:).Sample}';


SampleDir= [SampleInfo(1).ExperimentFolder '/' SampleInfo(1).AlignmentFolder ] ;

ainfo = load([SampleDir '/alignment_info']) ;
RefGenome = ainfo.ae.Genome ;

fr = fastaread(['/groups/kishony/Reference_Genomes/' RefGenome '/genome.fasta']) ;
GenomeLength=0;
ChrStarts=[];

ScafNames = {fr.Header} ;
for i=1:length(ScafNames)
    GenomeLength=GenomeLength+numel(fr(i).Sequence);
end

% get all coverage per bp and the mode coverage for each sample
[all_coverage_per_bp, coverage_modes] = get_all_coverage(SampleInfo, GenomeLength);
save('coveragematrix.mat', 'all_coverage_per_bp', 'coverage_modes', '-v7.3'); 
% 
% % normalize within sample by mode
% numisolates = size(all_coverage_per_bp,1);
% intervals = round(linspace(1, GenomeLength, 100)); 
% smoothwindow = intervals(2) - intervals(1) + 1; 
% mode_norm_cov = zeros(size(all_coverage_per_bp)); 
% 
% 
% fprintf('Normalizing within-sample by mode\n'); 
% for n = 1:numisolates
% 	fprintf('\nIsolate %i of %i\n', n, numisolates); 
% 	cur_isolate = double(all_coverage_per_bp(n,:)); 
% 	for i = 1:length(intervals) - 1
% 		if ~mod(i,100); fprintf('.'); end
% 		cur_interval = intervals(i):intervals(i+1); 
% 		mode_norm_cov(n,cur_interval) = smooth(cur_isolate(cur_interval),smoothwindow)./coverage_modes(n); 
% 	end
% 	clear cur_isolate; 
% end
% save('modenormalized.mat', 'smoothed_cov', '-v7.3'); 
% 
% 
% 
% fprintf('Normalizing within-sample by mode\n'); 
% for n = 1:numisolates
% 	fprintf('\nIsolate %i of %i\n', n, numisolates); 
% 	cur_isolate = double(all_coverage_per_bp(n,:)); 
% 	for i = 1:length(intervals) - 1
% 		if ~mod(i,100); fprintf('.'); end
% 		cur_interval = intervals(i):intervals(i+1); 
% 		mode_norm_cov(n,cur_interval) = smooth(cur_isolate(cur_interval),smoothwindow)./coverage_modes(n); 
% 	end
% 	clear cur_isolate; 
% end
% save('modenormalized.mat', 'mode_norm_cov', '-v7.3'); 
% 
% % normalize across samples 
% coverage_threshold = 10; 
% goodisolates = coverage_modes > coverage_threshold; 
% normalized_cov = zeros(numisolates, GenomeLength); 
% for b = 1:GenomeLength
% 	if ~mod(b,1000); fprintf('.'); end
% 	cur_median = median(mode_norm_cov(goodisolates,b)); 
% 	cur_norm = mode_norm_cov(:,b)./cur_median; 
% 	normalized_cov(:,b) = cur_norm; 
% end
% save('normalized.mat', 'normalized_cov', '-v7.3'); 
% 
% % bin for plotting 
