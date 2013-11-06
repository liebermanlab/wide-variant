% 2013 Hattie Chung
% This is written very poorly right now. 

% get type of mutation (N, S, I, P, etc.) 
types=[annotation_all.type];
typesmatrix=repmat(types',1,Nsample);

genes=[annotation_all.gene_num];

% genes with at least one mutation 
hasmutation_good = hasmutation(MutQual>qual_0,:); 
mut_gene_numbers = genes(sum(hasmutation_good,2)>0); 

% number of unique genes that get mutated (if length(y) < length(x), there
% are genes that get mutated more than once!) 
unique_genes = unique(mut_gene_numbers,'stable');
num_unique_genes=length(unique_genes);

genes_mut_cnt = zeros(num_unique_genes, 2); 

for g = 1:length(unique_genes)
    cnt = length(find(mut_gene_numbers==unique_genes(g)));
    genes_mut_cnt(g,:) = [unique_genes(g), cnt];
end

%% gene numbers that are mutated more than once

genenum_of_multiple_muts = unique_genes(genes_mut_cnt(:,2)>1); 
multimut_positions = cell(length(genenum_of_multiple_muts),1);
genotypes = cell(length(genenum_of_multiple_muts),1); 
genotype_degenerates = cell(length(genenum_of_multiple_muts),1); 

for m = 1:length(genenum_of_multiple_muts);
    this_gene = containers.Map(); 
    degenerates = []; 
    
    % get all structure entries that match gene number 
    mutpositions = find(genes==genenum_of_multiple_muts(m));
    multimut_positions{m} = mutpositions; 
    
    % ALL "genotypes" 
    mutcalls_m = Callsgood(mutpositions,:);
    
    for x = 1:size(mutcalls_m,2)
        g_x = mutcalls_m(:,x); 
        
        % add to key if no N's and not in current hash table 
        if ~isKey(this_gene,g_x) && isempty(strfind(g_x','N'))
            this_gene(g_x) = 1; 
        elseif ~isKey(this_gene,g_x)
            degenerates(end+1) = x;
        else
            this_gene(g_x) = this_gene(g_x)+1; 
        end
    end
    
    genotypes{m} = this_gene; 
    genotype_degenerates{m} = degenerates; 
end