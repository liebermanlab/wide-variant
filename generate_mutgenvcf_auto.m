function [MutGenVCF, calls, quals] =generate_mutgenvcf_auto(Positions,StrainDirs,SampleNames,ScafNames,qname,Parallel) 

global TEMPORARYFOLDER;



%Updated January 2012 to also return Call and Gen


MutGenVCF = VCF_struct ;

params = {} ;
for i=1:numel(StrainDirs)
    params{end+1} = {[StrainDirs{i} '/strain.vcf'],SampleNames{i}, ScafNames,Positions, TEMPORARYFOLDER} ;
end
run_parallel_matlab_commands('get_specific_line_multiple', params, qname, Parallel) ;


fprintf('Combining files...\n')
for i=1:numel(StrainDirs)
   % fprintf(1,'Sample: %g  ',i) ;
    vcf = load([TEMPORARYFOLDER '/vcfinfo_' SampleNames{i} '.mat']);
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
    delete([TEMPORARYFOLDER '/vcfinfo_' SampleNames{i} '.mat'])
end


[calls,gen] = build_call_including_indels(MutGenVCF) ;
quals = reshape(-1*[MutGenVCF.FQ],size(MutGenVCF)) ; 



return

