function [m, m_coding_strand] = div_mutation_type_probability_matrix(df, glength, ChrStarts, Sequences)

%Tami Lieberman May 2012
%Made from annotate_mutations_auto_gb and div_find_amino_acids

%df is coding sequences from annotate_mutations_auto_gb
%types='NSIPU';

Positions=p2chrpos([1:glength]',ChrStarts);

m=zeros(4,4,5); %original nt, new nt, type

m_coding_strand = zeros(4,4,2); %original nt on coding strand, new nt on coding stand, type

nts='atcg';
rc='tagc';

rv=[2 1 4 3];


an = zeros(size(Positions,1),1) ;
for i=1:numel(df)
    z = Positions(:,1)==i ;
    an(z) = genomic_position(df{i},Positions(z,2)) ;
    [~,Sequences{i}]=ismember(double(Sequences{i}),double(nts));
end




for i=1:size(Positions,1)

    ref=Sequences{Positions(i,1)}(Positions(i,2));
    
    
    if ref > 0
        
        nonref=1:4; nonref(ref)=[];
        
        if an(i) > 0 && an(i)==round(an(i)) % intragenic
            
            
            cdf = df{Positions(i,1)}(an(i)) ;
            
            if cdf.strand
                p = double(cdf.loc2) - double(Positions(i,2)) + 1;
                
            else
                p = double(Positions(i,2)) -double(cdf.loc1) + 1;
                if nts(ref) ~= cdf.Sequence(p)
                    disp(Positions(i))
                    disp(nts(ref))
                    disp(cdf.Sequence(p))
                end
            end
            
            
            aan = floor((p-1)/3) + 1 ;
            
            if numel(cdf.Sequence) >= aan*3
                
                ncn = p-(aan-1)*3 ;
                codon = cdf.Sequence(aan*3-2:aan*3) ;
                refAA= double(nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11));
                for j=nonref %for each nts
                    
                    
                    %find amino acid
                    if cdf.strand
                        codon(ncn)=rc(j);
                        coi=rv(ref); %coding strand original index
                        cni=rv(j);%coding strand new index
                    else
                        codon(ncn)=nts(j);
                        coi=ref;
                        cni=j;
                    end
                    
                    %is same?
                    if double(nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F', 'GENETICCODE', 11)) ~= refAA
                        m(ref,j,1) = m(ref,j,1)+1;
                        m_coding_strand(coi,cni,1) = m_coding_strand(coi,cni,1)+1;
                    else
                        m(ref,j,2) = m(ref,j,2)+1;
                        m_coding_strand(coi,cni,2) = m_coding_strand(coi,cni,2)+1;
                        
                    end
                    %setting ACGTonly to false prevents nt2aa from throwing an error for non-ACTG calls
                    %setting ALTERNATIVESTARTCODONS to false prevents nts from
                    %being called as alternative start codons
                end
                
            else
                m(ref,nonref,5)=m(ref,nonref,5)+1;
            end
        else
            
            p=0;
            if floor(an(i))>=1 % gene before position
                cdf = df{Positions(i,1)}(floor(an(i))) ;
                distance1 = Positions(i,2) - cdf.loc2;
                if cdf.strand==1
                    distance1 =distance1 * -1;
                end
                if distance1 <  0 & distance1 > -150
                    p=1;
                end
            end
            if ceil(an(i))<=length(df{Positions(i,1)}) % gene after position
                cdf = df{Positions(i,1)}(ceil(an(i))) ;
                distance2 = cdf.loc1 - Positions(i,2);
                if cdf.strand==0
                    distance2 = distance2 * -1;
                end
                if distance2 <  0 & distance2 > -150
                    p=1;
                end
            end
            
            
            if p==1
                m(ref,nonref,4)=m(ref,nonref,4)+1;
            else
                m(ref,nonref,3)=m(ref,nonref,3)+1;
            end
        end
    end
end


