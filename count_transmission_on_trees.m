function numtranmissions_lung = count_transmission_on_trees(genotypes, compartmenttimes, compartmentnumber, multiplestraininfection, mintimes)


%Each row in genotypes is a genotype, each column is a mutation
%Each row in compartmenttimes is a genotype, each column is a compartment,
%values are the number of times it was found in that compartment

%Tami Lieberman, July 2015


if nargin < 5
    mintimes=2;
end
source_compartment=1;



numtranmissions_lung=0;

genotypes_in_source=compartmenttimes(:,source_compartment)>=mintimes;
mutations_in_source=find(sum(genotypes(genotypes_in_source,:))>0);
    
genotypes_in_organ=compartmenttimes(:,compartmentnumber)>=mintimes;
mutations_in_organ=find(sum(genotypes(genotypes_in_organ,:),1)>0);

unaccounted_for_shared_mutations=mutations_in_organ(ismember(mutations_in_organ,mutations_in_source));

%disp(size(sum(genotypes(:,unaccounted_for_shared_mutations),2)>0))
%disp(size(genotypes_in_organ))
%disp(size(genotypes_in_source))


genotypes_containing_shared_mutation=sum(genotypes(:,unaccounted_for_shared_mutations),2)>0 & (genotypes_in_organ | genotypes_in_source);


if ~multiplestraininfection
    %deal with two edge cases of
    %(1) ancestor of patient found in organ
    if sum(sum(genotypes(genotypes_in_organ,:),2)==0)>0 %ancestor in organ
        numtranmissions_lung=numtranmissions_lung+1;
        %(2) where the ancestor isn't found in organ but a unique descendent of it is
    elseif  sum(genotypes_in_organ & ~genotypes_containing_shared_mutation) > 0
        numtranmissions_lung=numtranmissions_lung+1;
    end
end



while numel(unaccounted_for_shared_mutations)>0
    numtranmissions_lung=numtranmissions_lung+1;
    %find genotypes that contribute to these
    %unshared mutations
    %pick genotypes with fewest mutations first
    %first look for genotypes in organ 
    
    genotypes_containing_shared_mutation=find(sum(genotypes(:,unaccounted_for_shared_mutations),2)>0 & (genotypes_in_organ));
    if isempty(genotypes_containing_shared_mutation)
        genotypes_containing_shared_mutation=find(sum(genotypes(:,unaccounted_for_shared_mutations),2)>0 & (genotypes_in_source));
        fprintf(1,'blah');
    end
    minmuts=min(sum(genotypes(genotypes_containing_shared_mutation,:),2));
    mingenotype=genotypes_containing_shared_mutation(find(sum(genotypes(genotypes_containing_shared_mutation,:),2)==minmuts,1));
 
    unaccounted_for_shared_mutations(ismember(unaccounted_for_shared_mutations,find(genotypes(mingenotype,:))))=[];
    
    
end
