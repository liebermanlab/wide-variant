function [compartmentTimes, compartmentFreqs] = subsample_from_each_compartment(isolateCompartments,observedMatrix,genotypeMatrix)

%Tami Lieberman, August 2015
%This function does several things:
% 1) Removes locations with fewer than 5 samples
% 2) Subsamples each location to the minimum of remaining samples
% 3) Calculates compartmentTimes and compartmentFreqs on subsampled data

%%
global compartments

%% subsample to make a smaller observed matrix, based on isolate_compartments

%decide how many to take from each compartment
numPerCompartment=zeros(numel(compartments),1);
for c=1:numel(compartments)
    numPerCompartment(c)=sum(isolateCompartments==c);
end
goodCompartments=find(numPerCompartment>6);
toTakeFromEachCompartment=min(numPerCompartment(goodCompartments));

%fprintf(1,[num2str(toTakeFromEachCompartment) ' ']);

%randomly choose samples to take
newIsolateComparments=zeros(size(isolateCompartments));
for i=1:numel(goodCompartments)
    newIsolateComparments(datasample(find(isolateCompartments==goodCompartments(i)),toTakeFromEachCompartment,'Replace',false))=goodCompartments(i);
end

%make new observed matrix
observedMatrix=observedMatrix(:,newIsolateComparments>0);
newIsolateComparments=newIsolateComparments(newIsolateComparments>0);

%disp('lakdj')
%disp(newIsolateComparments)

%% assign to genotypes to make compartmentFreqs and compartmentTimes

mutsPerGenotype=sum(genotypeMatrix,2);
[~,genotypesByNumberMutations]=sort(mutsPerGenotype,'descend');

compartmentFreqs=zeros(length(genotypesByNumberMutations)+1,numel(compartments)); %last element is ancestral
compartmentTimes=zeros(length(genotypesByNumberMutations)+1,numel(compartments)); %last element is ancestral


%go each isolate, one at a time
for i=1:numel(newIsolateComparments)
    
    freqs=observedMatrix(:,i)';

    freqsfound=[];
    
    %try to assign to genotypes, starting with things with most mutations
    for j=1:numel(genotypesByNumberMutations)
        
        genotypen=genotypesByNumberMutations(j);
        g=find(genotypeMatrix(genotypen,:)>0); %which mutations are in this genotype
        
%        disp(min(freqs(g)))
        if min(freqs(g))>.15 %if all mutations in genotype found above 15% frequency
            
            genotypefreq=min(freqs(g));
            freqs(g)=freqs(g)-genotypefreq; %subtract this genotype from observation

            %record for data
          %  disp(newIsolateComparments(i))
            compartmentFreqs(j,newIsolateComparments(i))=compartmentFreqs(j,newIsolateComparments(i))+genotypefreq;
            compartmentTimes(j,newIsolateComparments(i))=compartmentTimes(j,newIsolateComparments(i))+1;
            freqsfound(end+1)=floor(genotypefreq*100);
        end
    end
    
    %ancestral case
    if sum(freqs) < .8 & sum(freqsfound)< 85 %if there isn't a lot left
        genotypefreq=(100-sum(freqsfound))/100;
        freqsfound(end+1)=floor(genotypefreq*100);
        compartmentFreqs(end,newIsolateComparments(i))=compartmentFreqs(end,newIsolateComparments(i))+genotypefreq;
        compartmentTimes(end,newIsolateComparments(i))=compartmentTimes(end,newIsolateComparments(i))+1;
    end
    
end
