function generate_parsimony_tree(calls, names)


if ~exist('dnapars','file')
    !cp /Volumes/sysbio/KISHONY\ LAB/illumina_pipeline/Tami/scripts/phylip-3.69/exe/dnapars dnapars
end
    

timestamp=datestr(now, 'yyyy-mm-dd-HH-MM-SS');

%write input file
generate_phylip_input(calls, names, [timestamp '_infile.txt'])

%write parameter file
fid = fopen([timestamp '_optionfile.txt'],'w');
fprintf(fid, [timestamp '_infile.txt\n']);
fprintf(fid, 'f\n');
fprintf(fid, [timestamp '_out.txt\n']);
fprintf(fid, 'y\n');
fprintf(fid, 'f\n');
fprintf(fid, [timestamp '_out.tree\n\n']);
fclose(fid);

%run
fprintf(1, 'Generating a tree according to maximum parsimony...')
fid = fopen('temp.sh','w');
fprintf(fid, ['! ./dnapars < ' timestamp '_optionfile.txt > ' timestamp '_outfile.txt\n']);
fclose(fid);
! chmod +x temp.sh
! ./temp.sh

