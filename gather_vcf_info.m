function [allcalls, allquals] =gather_vcf_info(Positions,StrainDirs,SampleNames,ScafNames,ChrStarts, qname,Parallel) 

global TEMPORARYFOLDER;



%Updated January 2012 to also return Call and Gen


allquals=zeros(size(Positions,1),numel(StrainDirs));
allcalls=char(allquals);

params = {} ;
for i=1:numel(StrainDirs)
    params{end+1} = {[StrainDirs{i} '/strain.vcf'],SampleNames{i}, ScafNames,ChrStarts, Positions, TEMPORARYFOLDER} ;
end
run_parallel_matlab_commands('gather_vcf_info_single_sample', params, qname, Parallel) ;





fprintf('Combining files...\n')
for i=1:numel(StrainDirs)
    fprintf(1,'Sample: %g  ',i) ;
    vcf = load([TEMPORARYFOLDER '/vcfinfo_' SampleNames{i} '.mat']);
    allcalls(:,i)=vcf.calls;
    allquals(:,i)=vcf.quals;
end

allquals=allquals*-1;

%delete things only after everything is loaded, helps for troubleshooting
for i=1:numel(StrainDirs)
    delete([TEMPORARYFOLDER '/vcfinfo_' SampleNames{i} '.mat'])
end


return

