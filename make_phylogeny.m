outgroup = 1; 
mainfolder = '../../';
RefGenomeNameUnedited = RefGenome;
% RefGenomeNameUnedited = 'Paeruginosa_DK2'; 
% TreeSampleNames = SampleNames; 
all_positions = 1:size(Calls,1); 

quality_positions = MutQual>qual_0 ;
quality_positions_chromosomal = p(quality_positions); 
%calls_for_tree = Calls(quality_positions,:); 
%dontplotisolates=[99 100];
toplot=ones(1,numel(names));
%toplot(dontplotisolates)=0;

%% fill in calls


calls_for_tree = Calls;
calls_for_tree=calls_for_tree(quality_positions,(~mixedisolates & toplot));


if outgroup == 1
    % ADD REFERENCE AT THESE POSITIONS
    [outgroup_nts, chromosome_name] = extract_outgroup_mutation_positions(mainfolder, RefGenomeNameUnedited, quality_positions_chromosomal); 
    calls_for_tree = [outgroup_nts', calls_for_tree]; 
    TreeSampleNames{1} = 'Outgroup'; 
    TreeSampleNames(2:sum(~mixedisolates & toplot)+1) = {SampleInfo(~mixedisolates & toplot).Sample};
else
    TreeSampleNames = SampleNames(SampleNames,:); 
end

generate_parsimony_tree(calls_for_tree, TreeSampleNames); 

fprintf('\nDone with tree\n'); 

%% write [chromosome name, quality_positions_chromosomal, calls for all isolates 

treecounting_csv_fn = 'test_treecounting/treecounting_chart.csv';
fid = fopen(treecounting_csv_fn, 'w'); 
num_muts = size(calls_for_tree,1); 
num_samples = length(TreeSampleNames); 

% write header
fprintf(fid, 'chr,pos');
for s = 1:num_samples
    fprintf(fid, ',%s', TreeSampleNames{s});
end
fprintf(fid,'\n');

% write data 
for m = 1:num_muts
    fprintf(fid, '%s,%i', chromosome_name, quality_positions_chromosomal(m)); 
    for s = 1:num_samples
        fprintf(fid, ',%s', calls_for_tree(m,s)); 
    end
    fprintf(fid, '\n'); 
end

fclose(fid);
fprintf('\nDONE!\n');