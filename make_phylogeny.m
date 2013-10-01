outgroup = 1; 
mainfolder = '/Volumes/sysbio/KISHONY LAB/illumina_pipeline';
% RefGenomeNameUnedited = 'Smaltophilia_K279';
RefGenomeNameUnedited = 'Paeruginosa_DK2'; 
TreeSampleNames = SampleNames; 
all_positions = 1:size(Calls,1); 
quality_positions = all_positions(MutQual>qual_0);
quality_positions_chromosomal = p(quality_positions); 
calls_for_tree = Calls(quality_positions,:); 

if outgroup == 1
    % ADD REFERENCE AT THESE POSITIONS
    outgroup_nts = extract_outgroup_mutation_positions(mainfolder, RefGenomeNameUnedited, quality_positions_chromosomal); 
    calls_for_tree = [calls_for_tree, outgroup_nts']; 
    TreeSampleNames{end+1} = 'Outgroup'; 
end

% generate_parsimony_tree(calls_for_tree, TreeSampleNames); 

% generate_parsimony_tree(NTs(maNT(sum(mutAF>0,2)>0,:)), SampleNames); 
% %the generated [timestamp]_out.tree file is best viewed in FigTree
% fprintf('\nDone with tree\n'); 
