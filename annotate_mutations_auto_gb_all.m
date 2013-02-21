function [an, df, mut] = annotate_mutations_auto_gb_all(Positions,Scaf,RefGenome)

tic; fprintf(1,'annotate_mutations... ') ;
an = zeros(size(Positions,1),1) ;
for i=1:length(Scaf)
    a=Scaf{i};f=find(a=='|',2,'last');fn = a(f(1)+1:f(2)-1) ; 
    fr = genbankread(['../Reference_Genomes/' RefGenome '/' fn '.gb']) ;
    df{i} = parse_all_locations_gb(fr.CDS, fr.Sequence) ;    
    z = Positions(:,1)==i ;
    an(z) = genomic_position(df{i},Positions(z,2)) ;
end

mut=[] ;
for i=1:size(Positions,1)
    if ~mod(i,100), fprintf(1,'.') ; end
    mut(i).gene_num = an(i) ;
    mut(i).scafold = Positions(i,1) ;
    mut(i).pos = Positions(i,2) ;

    nscf = Positions(i,1) ;
    if an(i)==round(an(i)) % intragenic
        cdf = df{nscf}(an(i)) ;
        mut(i).gene       = cdf.gene ;
        mut(i).protein    = cdf.product ;
        mut(i).protein_id = cdf.protein_id ;
        mut(i).strand     = cdf.strand ;
        mut(i).loc1       = cdf.loc1 ;
        mut(i).loc2       = cdf.loc2 ;
        mut(i).Sequence   = cdf.Sequence ;        
        if cdf.strand
            p =  Positions(i,2)-cdf.loc1 + 1;
        else
            p = -Positions(i,2)+cdf.loc2 + 1;
        end
        mut(i).nc_pos = p ;
        aan = floor((p-1)/3) + 1 ;
        ncn = p-(aan-1)*3 ;
        seq = cdf.Sequence(aan*3-2:aan*3) ;
        aa = [] ;
        for j=1:size(Call,2)
            if Call(i,j)~='N'
                seq(ncn) = Call(i,j) ;
                aa(end+1) = nt2aa(seq, 'ACGTOnly', false) ;
            end
        end
        mut(i).NonSyn = length(unique(aa))>1 ;
    else % intergenic
        if floor(an(i))>=1 % gene after
            cdf = df{nscf}(floor(an(i))) ;
            mut(i).fasta1.gene = cdf.gene ;
            % mut(i).fasta1.downstream = cdf.strand==0
        end
        if ceil(an(i))<=length(df{Positions(i,1)}) % gene before
            cdf = df{nscf}(ceil(an(i))) ;
            mut(i).fasta2.gene = cdf.gene ;
            % mut(i).fasta2.downstream = cdf.strand==1
        end
        mut(i).NonSyn = nan ;
    end
end
fprintf(1,'%6.1f min\n',toc/60) ;

return
