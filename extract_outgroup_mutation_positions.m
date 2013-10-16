function [ref_tree_positions, chr_name] = extract_outgroup_mutation_positions(ref_folder, RefGenome, chromosomal_positions)
    fastafile = [ref_folder '/Reference_Genomes/' RefGenome '/genome.fasta']; 
    
    fprintf('\nGetting outgroup sequences from reference genome\n\t%s\n', fastafile);
    fr = fastaread(fastafile) ;
    fasta_seq = fr.Sequence; 
    ref_tree_positions = fasta_seq(chromosomal_positions); 
    
    % HACK (HC 2013/10/04)
    divider_positions = strfind(fr.Header,'|'); 
    chr_name = fr.Header(divider_positions(3)+1:divider_positions(4)-1); 
end