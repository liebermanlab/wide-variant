function [an, df, mut, sequences] = annotate_mutations_auto_gb(Positions,Scaf,RefGenome)

nts='atcg';
rc='tagc';

tic; fprintf(1,'annotate_mutations...\n ') ;
an = zeros(size(Positions,1),1) ;
sequences={};
for i=1:length(Scaf)
    
    a=Scaf{i};f=find(a=='|',2,'last');fn = a(f(1)+1:f(2)-1) ;
    fr = genbankread(['../Reference_Genomes/' RefGenome '/' fn '.gb']) ;
    %SK: if no Sequence in gb file, use fasta
    if isempty(fr.Sequence)
        fr2=fastaread(['../Reference_Genomes/' RefGenome '/genome.fasta']);
        fr.Sequence=lower(fr2.Sequence);
    end
    
   
    genes = locustag_from_text(fr.CDS) ;
    genes = div_add_in_nonCDS(genes, fr.Features);
    df{i} = parse_all_locations_gb(genes, fr.Sequence) ;  %also reverses strands in this, sorts by position
    %sort by position on genome
    [~,sortedpositions]=sort([df{i}.loc1]);df{i}=df{i}(sortedpositions);
    
    
    z = Positions(:,1)==i ;
    an(z) = genomic_position(df{i},Positions(z,2)) ;
    
    sequences{end+1}=fr.Sequence;
    
end

%
mut=[] ;
for i=1:size(Positions,1)
    if ~mod(i,100), fprintf(1,'.') ; end
    mut(i).gene_num = an(i) ;
    mut(i).scafold = Positions(i,1) ;
    mut(i).pos = Positions(i,2) ;
    mut(i).ref=char(sequences{mut(i).scafold}(mut(i).pos)-32);
    
    nscf = Positions(i,1) ;
    if an(i)==round(an(i)) % intragenic
        cdf = df{nscf}(an(i)) ;
        
        %if product is more than one line, reshape
        if size(cdf.product,1)>1
            cdf.product=reshape(cdf.product',size(cdf.product,1)*size(cdf.product,2),1)';
        end
        mut(i).gene       = cdf.gene ;
        mut(i).protein    = cdf.product ;
        mut(i).protein_id = cdf.protein_id ;
        mut(i).strand     = cdf.strand ; %1 indicates reverse strand, 0 forward strand
        mut(i).loc1       = cdf.loc1 ;
        mut(i).loc2       = cdf.loc2 ;
        mut(i).Sequence   = cdf.Sequence ;
        mut(i).note   = cdf.note ;
        mut(i).locustag   = cdf.locustag ;
        mut(i).translation = cdf.translation;
                
        
        %fix frame annotation problems specific to Bdolosa genome
        if strcmp(RefGenome,'Bdolosa')
           % disp(mut(i).locustag)
            frame=ana_Bdolosa_checkframe(mut(i).locustag, Positions(i,2), double(mut(i).loc1));
        end
        
        mut(i).translation=nt2aa(mut(i).Sequence, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11) ;

        if frame~=0
            mut(i).Sequence=mut(i).Sequence(frame:end);
            if mut(i).strand==0
                mut(i).loc1 = cdf.loc1 + (frame - 1);
            else
                mut(i).loc2 = cdf.loc2 + (frame - 1);
            end
            mut(i).translation=nt2aa(mut(i).Sequence, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11) ;

        end
        
            
        if mut(i).strand
            p = double(mut(i).loc2) - double(Positions(i,2)) + 1;
        else
            p = double(Positions(i,2)) -double(mut(i).loc1) + 1;
        end
        mut(i).nt_pos = p ;
        
        aan = floor((p-1)/3) + 1 ;
        ncn = p-(aan-1)*3 ;
        codons=cell(4,1);
        AA='';
        mut(i).aa_pos=aan;
        
        if numel(mut(i).Sequence) >= aan*3 & mut(i).translation > 1;
            codon = mut(i).Sequence(aan*3-2:aan*3) ;
            for j=1:4 %for each nts
                if mut(i).strand
                    codon(ncn)=rc(j);
                else
                    codon(ncn)=nts(j);
                end
                codons{j}=codon;
                AA(j) = nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11) ;
                %setting ACGTonly to false prevents nt2aa from throwing an error for non-ACTG calls
                %setting ALTERNATIVESTARTCODONS to false prevents nts from
                %being called as alternative start codons
            end
        end
        
        
        mut(i).codons   = codons ;
        mut(i).AA   = AA ;
        mut(i).NonSyn = length(unique(AA))>1 ;
    else %intergenic
        if floor(an(i))>=1 % gene before position
            cdf = df{nscf}(floor(an(i))) ;
            mut(i).gene1 = cdf.gene ;
            mut(i).locustag1 = cdf.locustag ;
            mut(i).distance1 = Positions(i,2) - cdf.loc2;
            mut(i).protein1    = cdf.product ;
            if cdf.strand==1
                mut(i).distance1 = mut(i).distance1 * -1;
            end
        end
        if ceil(an(i))<=length(df{Positions(i,1)}) % gene after position
            cdf = df{nscf}(ceil(an(i))) ;
            mut(i).gene2 = cdf.gene ;
            mut(i).locustag2 = cdf.locustag;
            mut(i).distance2 = cdf.loc1 - Positions(i,2);
            mut(i).protein2    = cdf.product ;
            if cdf.strand==0
                mut(i).distance2 = mut(i).distance2 * -1;
            end
        end
        mut(i).NonSyn = nan;
    end
    
end
fprintf(1,'%6.1f min\n',toc/60) ;

return
