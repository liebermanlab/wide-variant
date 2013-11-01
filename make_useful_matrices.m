% Create useful matrices and vectors

global ANALYZE_DIVERSITY
[~, maNT, minorNT] = div_major_allele_freq(counts);

% ancnti = ancestral nucleotide in 0/1/2/3/4
% ancnt = ancestral nucleotide in N/A/T/C/G

if referenceisancestor > 0
    [ancnt, chromosome_name] = extract_outgroup_mutation_positions('/Volumes/sysbio/KISHONY LAB/illumina_pipeline', RefGenome, p); 
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

% for deep sequencing
if ANALYZE_DIVERSITY
    diversemutation=div_test_thresholds(counts,strict_parameters, coveragethresholds, CONTROLSAMPLE);
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
