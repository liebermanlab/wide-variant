%% Get StrainNames and ScafNames 
%***

vcf_q = 'all_2h' ;
Parallel = true ;

StrainNames = read_strain_names ;
NStrain = length(StrainNames) ;

RefGenome = {} ;
StrainDirs = {} ;

for i=1:NStrain
    StrainDirs{i} = ['../' StrainNames(i).ExperimentFolder '/' StrainNames(i).Sample '/' StrainNames(i).AlignmentFolder ] ;
    ainfo = load([StrainDirs{i} '/alignment_info']) ;
    RefGenome{i} = ainfo.ae.Genome ;
end

RefGenome = unique(RefGenome) ;
if length(RefGenome)>1
    error('Must compare samples aligned to the same reference genome')
end

RefGenome = RefGenome{1} ;

fr = fastaread(['../Reference_Genomes/' RefGenome '/genome.fasta']) ;

ScafNames = {fr.Header} ; 

for i=1:length(ScafNames) 
    f=find(ScafNames{i}==' ',1) ;
    ScafNames{i} = ScafNames{i}(1:f-1) ;
end

%% Build position list
switch 2 %***
    case 1
        [Positions,Genotypes] = generate_positions_genotypes_auto(StrainDirs, ScafNames, 50000) ;
        save positions_genotypes Positions Genotypes StrainNames ScafNames
    case 2
        load positions_genotypes
end


%% build VCF
switch 2 %***
    case 1
        MutGenVCF = generate_mutgenvcf_auto(Positions,StrainDirs,ScafNames,vcf_q,Parallel) ;                
        save mutgenvcf MutGenVCF
    case 2
        load mutgenvcf
end


%% Build Call

[Call,Gen] = build_call(MutGenVCF) ;

%% Remove identical lines
Call2 = Call ; 
%***
Call2(Call=='n' | Call=='D' | Call=='I') = 'N' ;
skip_v = zeros(size(MutGenVCF,1),1);
for k=1:size(Call2,1)
    skip_v(k) = length(unique([Call2(k,:), 'N']))<=2 ;
end
Positions3 = Positions(~skip_v,:) ;
MutGenVCF3 = MutGenVCF(~skip_v,:) ;
Call3 = Call2(~skip_v,:) ;
fprintf('Remove %g mutations with identical, or no calls, or deletions, or insertions\n',sum(skip_v)) ;


%% Annotate mutations
[Fasta3, df, Mut3] = annotate_mutations_auto_gb(Positions3,Call3,ScafNames,RefGenome) ;   
Qual3 = reshape([MutGenVCF3.qual],size(MutGenVCF3)) ;

save mutgenvcf_3 MutGenVCF3 Positions3 Fasta3 Mut3 Call3 Qual3 StrainNames NStrain ScafNames RefGenome

