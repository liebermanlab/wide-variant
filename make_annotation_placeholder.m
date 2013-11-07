function placeholder_annot = make_annotation_placeholder(node_name, contig_seq)
    loc1 = 0; 
    loc2 = length(contig_seq); 
    loc_str = sprintf('%i..%i', loc1, loc2); 
    placeholder_annot = struct('location', loc_str, ... % make it length of contig 
                                'gene', [], ...
                                'product', node_name,...
                                'codon_start', '1', ...
                                'indices', [loc1 loc2], ...
                                'protein_id', '', ...
                                'db_xref', '', ...
                                'note', '', ...
                                'translation', [], ...
                                'text', '', ...
                                'locustag', '', ...
                                'loc1', [loc1], ...
                                'loc2', [loc2], ...
                                'strand', 1, ...
                                'Sequence', contig_seq); 
end