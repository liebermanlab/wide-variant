function [u, l, expected_trials_matrix, percentN_matrix, original_matrix, expected, expectedI, expectedP, cds_p_types] = div_permute_expected(num_trials, actual_mutations,p_types, cds_p_types)


% p_types is of size 4 x 4 x 5  (N/S/I/P/U)

percentN_matrix=p_types(:,:,1)./(p_types(:,:,1)+p_types(:,:,2));

% find nt freq in genome and in cds

nt_freq=zeros(4,1);
cds_nt_freq=zeros(4,1);
for i=1:4
    nt_freq(i)=sum(p_types(i,:));
    cds_nt_freq(i)=sum(p_types(i,:,1))+sum(p_types(i,:,2));
end
nt_freq=nt_freq./sum(nt_freq(:));
cds_nt_freq=cds_nt_freq./sum(cds_nt_freq(:));


% generate numtrials mutation incidence matrices, picking with replacement
% bootstrapped_matrix is of dimension num_trials x 4 x 4
% bootstrapped_matrix is normalized by total mutations and original nt freq
[bootstrapped_matrix, original_matrix] = div_expected_matrix_permutation(num_trials, actual_mutations, nt_freq);


normalized_matrix=zeros(4,4);
% normalize bootstrapped_matrix to cds_nt_freq
for i=1:4
    bootstrapped_matrix(:,i,:)=bootstrapped_matrix(:,i,:).*cds_nt_freq(i);
    cds_normalized_matrix(i,:)=original_matrix(i,:).*cds_nt_freq(i);
    normalized_matrix(i,:)=original_matrix(i,:).*nt_freq(i);
end



%combine each bootstrapped matrix with NS ratio to get a simple expected
%percentN for each matrix

percentN=zeros(num_trials,1);
for i=1:num_trials
    
    %first normalize matrices so that the sum is 1
    occurence_matrix=percentN_matrix.*squeeze(bootstrapped_matrix(i,:,:)./sum(bootstrapped_matrix(i,:)));
    occurence_matrix(isnan(occurence_matrix))=0;
    %then combine
    percentN(i)=sum(occurence_matrix(:));

end

expected_occurence=percentN_matrix.*cds_normalized_matrix./sum(cds_normalized_matrix(:));
expected_occurence(isnan(expected_occurence))=0;
expected=sum(expected_occurence(:))/(1-sum(expected_occurence(:)));


percentI_matrix=(p_types(:,:,3)+ p_types(:,:,4))./(sum(p_types,3)); 
Ioccurence=percentI_matrix.*normalized_matrix./sum(normalized_matrix(:)); Ioccurence(isnan(Ioccurence))=0;
percentP_matrix=p_types(:,:,4)./(sum(p_types,3)); percentP_matrix(isnan(percentP_matrix))=0;
Poccurence=percentP_matrix.*normalized_matrix./sum(normalized_matrix(:)); Poccurence(isnan(Poccurence))=0;

expectedI=sum(Ioccurence(:))/(1-sum(Ioccurence(:)));
expectedP=sum(Poccurence(:))/(1-sum(Poccurence(:)));


u=zeros(1,1000);
l=zeros(1,1000);

bins=1000;


%run simulation num_trials time for each bin size
%each set of simulations for a different number of mutations is run with the same num_trials matrices



expected_trials_matrix=zeros(num_trials,bins);

for n=1:bins
    
    %calculate N and S for each distribution
    
    N=nan(num_trials,1);
    
    for t=1:num_trials
       r=rand(n,1);
       N(t)=sum(r<=percentN(t));
    end
    N=sort(N);
    u(n)=N(floor(.975*num_trials)+1)/(n-N(floor(.975*num_trials+1))); %N/S
    l(n)=N(floor(.025*num_trials))/(n-N(floor(.025*num_trials))); %N/S
    expected_trials_matrix(:,n)=N./(n-N);
end





end