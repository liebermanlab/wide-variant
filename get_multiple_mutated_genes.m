% get type of mutation (N, S, I, P, etc.) 
types=[annotation_all.type];
typesmatrix=repmat(types',1,Nsample);

genes=[annotation_all.gene_num];

% num_genes_with_at_least_one_mut
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
multimut_gene_num = unique_genes(genes_mut_cnt(:,2)>1); 
for m = 1:length(multimut_gene_num)
    temp = find(genes==multimut_gene_num(m)); 
    multiple_gene_index(m) = temp(1); 
    annotation_all(temp(1))
end