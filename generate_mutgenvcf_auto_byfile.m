function MutGenVCF =generate_mutgenvcf_auto_byfile(Positions,StrainFiles,ScafNames,qname,Parallel) 

%generate_mutgenvcf_auto(target_dir, Positions,SNPFileName,ScafNames,RefGenome)
MutGenVCF = VCF_struct ;

params = {} ;
for i=1:length(StrainFiles)
    params{i} = {StrainFiles{i},ScafNames,Positions,sprintf('tmp%g',i)} ;
end

run_parallel_matlab_commands('get_specific_line', params, qname, Parallel) ;

for i=1:length(StrainFiles)
    fprintf(1,'Strain: %g  ',i) ;
    vcf = load(sprintf('tmp%g',i)) ;
    for k=1:size(Positions,1) ;
        j = find(vcf.indx==k) ;
        if isempty(j)
            MutGenVCF(k,i).qual = nan ;
            MutGenVCF(k,i).pos = nan ;
            MutGenVCF(k,i).gen2 = nan ;
            MutGenVCF(k,i).FQ = nan ;
        else
            MutGenVCF(k,i) = read_vcf_line(vcf.lins{j}) ;
        end
    end
end

return

