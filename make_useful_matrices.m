% Create useful matrices and vectors

[~, maNT, minorNT] = div_major_allele_freq(counts);

% ancnti = ancestral nucleotide in 0/1/2/3/4
% ancnt = ancestral nucleotide in N/A/T/C/G

if ancestoriscontrol>0
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

% for deep sequencing
diversemutation=div_test_thresholds(counts,strict_parameters, coveragethresholds, CONTROLSAMPLE);
% for isolates
fixedmutation=((Calls~=ancnt_m) & (ismember(Calls,acceptabletypes)) & repmat(MutQual,1,Nsample)>=qual_0);
% either type of mutation (fixed or diverse) 
hasmutation= fixedmutation | diversemutation; %has mutation if diverse or call(from vcf file) ~= anct

%Mutant allele frequency
[~, mutantNT]=mutant_frequency(counts, hasmutation, ancnti, Calls);

mut_freq = get_fixed_mut_freq(Calls, MutQual, qual_0, ancnti);

Callsgood = Calls(MutQual>qual_0, :); 
countsgood = counts(:,MutQual>qual_0,:); 
isN = get_ambiguous_calls(Callsgood);