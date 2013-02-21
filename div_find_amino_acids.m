function [refnts, aminoacids, pintergenic, ppromoter] = div_find_amino_acids(Positions, df, promoterdistance)

%created from annotate_mutations_auto_gb
%designed to be called from div_calculate_expected_matrix
%Tami Lieberman May 2012

nts='atcg';
rc='tagc';

an = zeros(size(Positions,1),1) ;
for i=1:numel(df)
    z = Positions(:,1)==i ;
    an(z) = genomic_position(df{i},Positions(z,2)) ;
end

refnts = [] ;
aminoacids = [] ;
pintergenic=0;
ppromoter=0;



for i=1:size(Positions,1)
    
    if ~mod(i,1000), fprintf(1,'.') ; end
    
    if an(i) > 0 && an(i)==round(an(i)) % intragenic
        
        
        cdf = df{Positions(i,1)}(an(i)) ;
        
        if cdf.strand
            p = double(cdf.loc2) - double(Positions(i,2)) + 1;
        else
            p = double(Positions(i,2)) -double(cdf.loc1) + 1;
        end
        aan = floor((p-1)/3) + 1 ;
        
        if numel(cdf.Sequence) >= aan*3 & find(double(nts)==double(cdf.Sequence(p)));
            ncn = p-(aan-1)*3 ;
            AAs = [ 0 0 0 0 ];
            codon = cdf.Sequence(aan*3-2:aan*3) ;
            for j=1:4 %for each nts
                if cdf.strand
                    codon(ncn)=rc(j);
                else
                    codon(ncn)=nts(j);
                end
                AAs(j) = double(nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F')) ;
                %setting ACGTonly to false prevents nt2aa from throwing an error for non-ACTG calls
                %setting ALTERNATIVESTARTCODONS to false prevents nts from
                %being called as alternative start codons
            end
            if cdf.strand
                refnts = [refnts find(double(rc)==double(cdf.Sequence(p)))] ;
            else
                refnts = [refnts find(double(nts)==double(cdf.Sequence(p)))] ;
            end
            aminoacids = [ aminoacids ; AAs ] ;
        end
    else
        pintergenic=pintergenic+1;
        
        if floor(an(i))>=1 % gene before position
            cdf = df{Positions(i,1)}(floor(an(i))) ;
            distance1 = Positions(i,2) - cdf.loc2;
            if cdf.strand==1
                distance1 =distance1 * -1;
            end
        end
        if ceil(an(i))<=length(df{Positions(i,1)}) % gene after position
            cdf = df{Positions(i,1)}(ceil(an(i))) ;
            distance2 = cdf.loc1 - Positions(i,2);
            if cdf.strand==0
                distance2 = distance2 * -1;
            end  
        end        
        if (distance1 <  0 & distance1 > -150) | (distance2 <  0 & distance2 > -150)
            ppromoter=ppromoter+1;
        end
    end
    
end

pintergenic=pintergenic/size(Positions,1);
ppromoter=ppromoter/size(Positions,1);


return
