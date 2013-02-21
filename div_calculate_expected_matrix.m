function [expectedN, expectedS, pmatrix, expected_intergenic, expected_promoter, ntfreq] = div_calculate_expected_matrix(Npositions, cds, actual_mutations, GLength, ChrStarts, promotersize) % Ref, Scaf, GLength, ChrStarts)

% May 2012, Tami Lieberman

nts=double('atcgr');
rc=double('tagc');
ntfreq=zeros(4,1);


fmatrix=zeros(4,4);
pmatrix=zeros(4,4);

for i=1:numel(actual_mutations)
    %if in a gene
    if floor(actual_mutations(i).gene_num)==actual_mutations(i).gene_num
        
        %find old nt
        if actual_mutations(i).strand
            rc_original = double(actual_mutations(i).Sequence(actual_mutations(i).nt_pos));
            original_index = find(rc==rc_original,1);
        else
            original= double(actual_mutations(i).Sequence(actual_mutations(i).nt_pos));
            original_index = find(nts==original,1);
        end
        
        %find new nt
        pair=actual_mutations(i).nts+32; %move to lowercase to match sequence
        pair(pair==nts(original_index))=[]; %remove original
        for j=1:numel(pair) %if more than one nts count each occurence as if it mutated from the reference
            new_index=find(nts==pair(j));
            fmatrix(new_index,original_index)=fmatrix(new_index,original_index)+1;
        end
    end
    
end

fprintf('Done creating matrix\n')

%find random positions on genomes
p = floor(rand(Npositions,1)*GLength)+1;
positions=p2chrpos(p,ChrStarts);
[refnts, aminoacids, expected_intergenic, expected_promoter] = div_find_amino_acids(positions,cds, promotersize);

%normalize occurence matrix to turn into probability matrix
%this requires considering frequency of each nt in reference genome
for i=1:4
    ntfreq(i)=sum(refnts==i);
    pmatrix(:,i)=fmatrix(:,i)/ntfreq(i);
end
pmatrix=pmatrix./sum(pmatrix(:));

fprintf('Done normalizing matrix\n')

%find new nts based on probability matrix
Npositions=size(aminoacids,1); %removes positions
x=rand(Npositions,1); %random number < 1
randoms=[x x x x]';

mutatematrix=pmatrix(:,refnts);
[sorted,sortednts]=sort(mutatematrix);
mutatematrix=cumsum(sorted);

n=sum(mutatematrix<=randoms)+1;
n(n==5)=1;
newnts=sortednts(sub2ind(size(sortednts), n, 1:Npositions));

mutates= find(newnts~=refnts);
 

%count types of mutations

fprintf(['Testing N/S at ' num2str(numel(mutates)) ' nt positions'])

N = 0 ; S = 0 ;

confirmmatrix=zeros(4,4);

for i=1:numel(mutates)
    index=mutates(i);
    confirmmatrix(newnts(index),refnts(index)) = confirmmatrix(newnts(index),refnts(index)) + 1;
    if aminoacids(index,refnts(index))==aminoacids(index,newnts(index))
        S=S+1;
    else
        N=N+1;
    end
            
end


%this matrix should be compared with frequency matrix, not probability matrix
confirmmatrix=confirmmatrix/sum(confirmmatrix(:));

expectedN=N;
expectedS=S;

end
