% Create useful matrices and vectors

global ANALYZE_DIVERSITY

% ancnti = ancestral nucleotide in 0/1/2/3/4
% ancnt = ancestral nucleotide in N/A/T/C/G

if referenceisancestor > 0
    ancnt = extract_outgroup_mutation_positions('~/Dropbox/', RefGenome, positions); 
    [~,ancnti]=ismember(ancnt',NTs);
    ancnt=ancnt';
elseif ancestoriscontrol>0
    %Ancestral nucleotide at each position
    disp(['Using sample ' num2str(ancestoriscontrol) ' of major alleles as ancestor...\n']);
    [ancnt, ancnti] = ancestorfromsample(Calls,CONTROLSAMPLE);
elseif ancestorismode==1
    disp('Using mode of major alleles as ancestor...\n');
    ancnti=mode(maNT,2);
    temp=ancnti; temp(ancnti>4 | ancnti<1)=5; NTsN='ATCGN';
    ancnt=NTsN(temp)';
else
    disp('Error: choose setting for ancestor')
    stop
end
    
% expand to match # isolates    
ancnt_m=repmat(ancnt,1,Nsample); 
ancnti_m=repmat(ancnti,1,Nsample);

% Hasmutation


% for isolates
fixedmutation=((Calls~=ancnt_m) & (ismember(Calls,acceptabletypes)) & repmat(MutQual,1,Nsample)>=qual_0);

load(['cov_' run_postfix]);
coveragethresholds=coveragethresholds(:,goodsamples);
COVERAGETHRESHOLDS=coveragethresholds;



% for deep sequencing
if ANALYZE_DIVERSITY
    diversemutation=div_test_thresholds(counts,STRICT_PARAMETERS, COVERAGETHRESHOLDS, CONTROLSAMPLE);
else
    diversemutation=zeros(size(fixedmutation));
    
end
% either type of mutation (fixed or diverse) 
hasmutation= fixedmutation | diversemutation; %has mutation if diverse or call(from vcf file) ~= anct

%Mutant allele frequency
[mutAF, mutantNT]=mutant_frequency(counts, hasmutation, ancnti, Calls);


mut_freq = get_fixed_mut_freq(Calls, MutQual, qual_0, ancnti);

Callsgood = Calls(MutQual>qual_0, :); 
countsgood = counts(:,MutQual>qual_0,:); 
isN = get_ambiguous_calls(Callsgood);




%Call indels

[x,y]=find(dindel_freqs>.8);
candidate_indel_p=unique(x);
%for i=1:numel(candidate_indel_p)
    