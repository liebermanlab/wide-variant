function mutations = annotate_mutations_metagenomics(positions)

%ideally we should add marker annotations within this function


nts='atcg';

fasta=fastaread('trimmed.fna');


mutations=[] ;

for i=1:size(Positions,1)
   
    mutations(i).marker = Positions(i,1) ;
    mutations(i).pos = Positions(i,2) ;
    
    %SEARCH
    
    mutations(i).name = fasta(mutations(i).marker).Header ;
    mutations(i).sequence = fasta(mutations(i).marker).Sequence;
    
    
    p = mutations(i).pos;
    
    aan = floor((p-1)/3) + 1 ;
    ncn = p-(aan-1)*3 ;
    codons=cell(4,1);
    AA='';
    codon = cdf.Sequence(aan*3-2:aan*3) ;
    
    for j=1:4 %for each nts
        codon(ncn)=nts(j);
        codons{j}=codon;
        AA(j) = nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F') ;
        %setting ACGTonly to false prevents nt2aa from throwing an error for non-ACTG calls
        %setting ALTERNATIVESTARTCODONS to false prevents nts from
        %being called as alternative start codons
    end
    
    mutations(i).codons   = codons ;
    mutations(i).AA   = AA ;

end

