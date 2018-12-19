function [allcalls, allquals] =gather_vcf_info(Positions,StrainDirs,SampleNames,ScafNames,ChrStarts, qname,Parallel)

global TEMPORARYFOLDER;



%Updated January 2012 to also return Call and Gen


allquals=zeros(size(Positions,1),numel(StrainDirs));
allcalls=char(allquals);



save('for_gather_vcf_info', 'ScafNames','ChrStarts', 'Positions', 'TEMPORARYFOLDER', '-v7.3')


cmds = {} ;
for i=1:numel(StrainDirs)
    
    
    if ~(exist([TEMPORARYFOLDER '/vcfinfo_' SampleNames{i} '.mat'],'file'))
        vcfpath=[StrainDirs{i} '/strain.vcf'];
        cmds{end+1} = ['matlab -r "path(' char(39) '/scratch/users/tami/illumina_pipeline_c3ddb/' char(39) ',path); gather_vcf_info_single_sample(' char(39) vcfpath char(39) ',' char(39) SampleNames{i} char(39) ');"'];
        
    end
    %     if numel(params) >800
    %         run_parallel_matlab_commands('gather_vcf_info_single_sample', params, qname, Parallel) ;
    %         params={};
    %     end
    
end
run_parallel_unix_commands_fast(cmds,qname,Parallel,{'.'});





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

