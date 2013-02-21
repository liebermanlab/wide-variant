%% One test

actual_mutations=combined_annotation_goodpos;
genomeNmatrix=percentN_matrix;


%Tami Lieberman June 2012

NTs='ATCG';
nts='atcg';

actual_matrix=zeros(4,4,2);


for i=1:numel(actual_mutations)
    
    if numel(actual_mutations(i).AA)==4
        
        original_index = find(NTs==actual_mutations(i).ref,1);
        
        %find new nt
        pair=actual_mutations(i).nts; %move to lowercase to match sequence
        pair(pair==actual_mutations(i).ref)=[]; %remove original
        
        refAA=actual_mutations(i).AA(original_index);
        
        for j=1:numel(pair) %if more than one nts count each occurence as if it mutated from the reference
            new_index=find(NTs==pair(j));
            newAA=actual_mutations(i).AA(new_index);
            if newAA ~= refAA
                actual_matrix(original_index,new_index,1)=actual_matrix(original_index,new_index,1)+1;
            else
                actual_matrix(original_index,new_index,2)=actual_matrix(original_index,new_index,2)+1;
            end
        end
    end
    
end

%compress into 6 elements

%AT
%AC
%AG
%GC
%GT
%GA

actual=zeros(6,2);
genomeN=zeros(6,1);

actual(1,:)=squeeze(actual_matrix(1,2,:)+actual_matrix(2,1,:));
actual(2,:)=squeeze(actual_matrix(1,3,:)+actual_matrix(2,4,:));
actual(3,:)=squeeze(actual_matrix(1,4,:)+actual_matrix(2,3,:));
actual(4,:)=squeeze(actual_matrix(4,3,:)+actual_matrix(3,4,:));
actual(5,:)=squeeze(actual_matrix(4,2,:)+actual_matrix(3,1,:));
actual(6,:)=squeeze(actual_matrix(4,1,:)+actual_matrix(3,2,:));

genomeN(1,:)=(genomeNmatrix(1,2)+genomeNmatrix(2,1))/2;
genomeN(2,:)=(genomeNmatrix(1,3)+genomeNmatrix(2,4))/2;
genomeN(3,:)=(genomeNmatrix(1,4)+genomeNmatrix(2,3))/2;
genomeN(4,:)=(genomeNmatrix(4,3)+genomeNmatrix(3,4))/2;
genomeN(5,:)=(genomeNmatrix(4,2)+genomeNmatrix(3,1))/2;
genomeN(6,:)=(genomeNmatrix(4,1)+genomeNmatrix(3,2))/2;

actualN=actual(:,1)./(actual(:,1)+actual(:,2));



sd=sqrt((genomeN.*(1-genomeN))./sum(actual,2));

lower=genomeN-2.*sd;
upper=genomeN+2.*sd;


%% Another thing (temp.m)


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


