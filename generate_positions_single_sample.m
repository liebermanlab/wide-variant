function Positions = generate_positions_single_sample(StrainDir, SampleName, ScafNames, maxFQ, onlysnps, L, tempfolder)
% unlike generate_mutations_genotypes, this functions creates a list of mutated positions (chromosome, position) rather than a list of mutaions (chromosome, position, ref, alt).
% changes are mainly in not reading reference and alelle columns.

% an earlier version created a structure called Genotypes, this has since
% been removed

if nargin < 5
    tempfolder = '';
end


K = 0 ;
Positions = zeros(L,2) ;
Mnum = zeros(L,1) ;
include=zeros(L,1);

%fprintf(1,'Strain: %g  ',i) ;
vcf = read_vcf_file([StrainDir '/variant.vcf']) ;
if isstruct(vcf)
    for j=1:length(vcf) ;
        ScfN = find(strcmp(vcf(j).scaf,ScafNames)) ;
        pos = vcf(j).pos ;
        num = ScfN*1e8+pos ;
        Mindx=find(Mnum(1:K)==num,1) ;
        if ((isempty(Mindx)) & (numel(vcf(j).alt)==1 | onlysnps==0))
            K=K+1 ;
            Positions(K,:) = [ScfN, pos] ;
            Mnum(K) = num ;
            t = read_vcf_info(vcf(j).info) ;
            if t.FQ < maxFQ
                include(K)=1;
            end
        end
    end
    Positions=Positions(include>0,:);
else
    Positions=Positions(1,:);
end

save([tempfolder '/vcf_' char(SampleName)], 'Positions')
