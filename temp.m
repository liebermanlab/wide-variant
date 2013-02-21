
%Tami Lieberman June 2012, Neutral model


m=zeros(4,4,2); %original nt, new nt, type

nts='atcg';

codonsN=nan(64,3);

codons=nan(64,3);
index=1;
for p1=1:4
    
    for p2=1:4
        
        for p3=1:4
            
            codonsN(index,:)=[p1 p2 p3];
            codons(index,:)=[nts(p1) nts(p2) nts(p3)];
            index= index+1;
            
        end
    end
end
codons=char(codons);


for i=1:size(codons,1)
    
    
    
    codon=codons(i,:);
    
    refAA= double(nt2aa(codon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F'));
    
    if refAA ~=42
        
        for ncn=1:3
            
            ref=codonsN(i,ncn);
            nonref=1:4; nonref(ref)=[];
            
            
            newcodon=codon;
            
            for j=nonref %for each nts
                
                newcodon(ncn)=nts(j);
                
                %is same?
                if double(nt2aa(newcodon, 'ACGTOnly', false, 'ALTERNATIVESTARTCODONS','F')) ~= refAA
                    m(ref,j,1) = m(ref,j,1)+1;
                else
                    m(ref,j,2) = m(ref,j,2)+1;
                    
                end
                %setting ACGTonly to false prevents nt2aa from throwing an error for non-ACTG calls
                %setting ALTERNATIVESTARTCODONS to false prevents nts from
                %being called as alternative start codons
            end
            
        end
    end
end

c=m(:,:,1)./(m(:,:,1) +m(:,:,2) );
