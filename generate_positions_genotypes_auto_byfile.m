function [Positions,Genotypes] = generate_positions_genotypes_auto_byfile(StrainFiles, L)      
% unlikea generate_mutations_genotypes, this functions creates a list of mutated positions (chromosome, position) rather than a list of mutaions (chromosome, position, ref, alt). changes are mainly in not reading reference and alelle columns.  
%

Genotypes = struct('Mindx',{},'Mfreq',{}) ;
K = 0 ; 
Positions = zeros(L,2) ;
Mnum = zeros(L,1) ;

for i=1:length(StrainFiles)    
    fprintf(1,'Strain: %g  ',i) ;    
    vcf = read_vcf_file(StrainFiles{i}) ;     
    Genotypes(i).Mindx = zeros(length(vcf),1) ;
    Genotypes(i).Mfreq = zeros(length(vcf),1) ;
    for j=1:length(vcf) ;
        if ~mod(j,1000), fprintf(1,'.');  end
        ScfN = find(strcmp(vcf(j).scaf,ScafNames)) ;
        pos = vcf(j).pos ;
        num = ScfN*1e8+pos ;
        Mindx=find(Mnum(1:K)==num,1) ;
        if isempty(Mindx)
            K=K+1 ;
            Positions(K,:) = [ScfN, pos] ;
            Mnum(K) = num ;
            Mindx = K ; 
        end
        Genotypes(i).Mindx(j) = Mindx ;
        Genotypes(i).Mfreq(j) = vcf(j).gen2 ;
    end
    fprintf(1,'\n',i) ;
end

Positions = Positions(1:K,:) ;

return
