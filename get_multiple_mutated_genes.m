% get type of mutation (N, S, I, P, etc.) 
types=[annotation_all.type];
typesmatrix=repmat(types',1,Nsample);

genes=[annotation_all.gene_num];

% num_genes_with_at_least_one_mut
x = genes(sum(hasmutation,2)>0); 

% number of unique genes that get mutated (if length(y) < length(x), there
% are genes that get mutated more than once!) 
y=unique(x);