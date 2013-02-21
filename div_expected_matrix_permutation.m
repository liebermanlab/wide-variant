function  [matrices, pmatrix] = div_expected_matrix_permutation(numtrials, actual_mutations, ntfreq) % Ref, Scaf, GLength, ChrStarts)

% May 2012, Tami Lieberman

NTs=double('ATCGR');


fmatrix=zeros(4,4);
pmatrix=zeros(4,4);

pairs=zeros(2,numel(actual_mutations));
index=1;

for i=1:numel(actual_mutations)
    
    
    original_index = find(NTs==actual_mutations(i).ref,1);
    
    %find new nt
    pair=actual_mutations(i).nts; %move to lowercase to match sequence
    pair(pair==actual_mutations(i).ref)=[]; %remove original
    
    for j=1:numel(pair) %if more than one nts count each occurence as if it mutated from the reference
        new_index=find(NTs==pair(j));
        fmatrix(original_index, new_index)=fmatrix(original_index, new_index)+1;
        pairs(:,index)=[original_index new_index];
        index= index+ 1;
    end
    
end

pairs=pairs(:,1:index-1);
fprintf('Done creating original matrix\n')

matrices=zeros(numtrials,4,4);
for i=1:numtrials
    r=floor(rand(index-1,1)*(index-1))+1;
    newpairs=pairs(:,r);
    for j=1:size(newpairs,2);
        matrices(i, newpairs(1,j), newpairs(2,j))= matrices(i, newpairs(1,j), newpairs(2,j))+1;
    end
end


%normalize occurence matrix to turn into probability matrix
%this requires considering frequency of each nt in reference genome
for i=1:4
    matrices(:,i,:)=matrices(:,i,:)./ntfreq(i);
    pmatrix(i,:)=fmatrix(i,:)./ntfreq(i);
end

%normalize to have sum be 1
pmatrix=pmatrix./sum(pmatrix(:));
for i=1:numtrials
    matrices(i,:,:)=matrices(i,:,:)./sum(matrices(i,:));
end


end
