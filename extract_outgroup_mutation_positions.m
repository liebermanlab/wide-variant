function ref_tree_positions = extract_outgroup_mutation_positions(ref_folder, RefGenome, chromosomal_positions)
    fastafile = [ref_folder '/Reference_Genomes/' RefGenome '/genome.fasta']; 
    
    fprintf('\nGetting outgroup sequences from reference genome\n\t%s\n', fastafile);
    fr = fastaread(fastafile) ;
    fasta_seq = fr.Sequence; 
    ref_tree_positions = fasta_seq(chromosomal_positions); 
end