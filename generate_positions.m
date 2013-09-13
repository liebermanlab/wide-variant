function p = generate_positions(SampleDirs, SampleNames, GenomeLength, ScafNames, ChrStarts, maxFQ, onlysnps, L, parallel, jobsubmitoptions)      
% unlike generate_mutations_genotypes, this functions creates a list of mutated positions (chromosome, position) rather than a list of mutaions (chromosome, position, ref, alt).
% changes are mainly in not reading reference and alelle columns.  

% an earlier version created a structure called Genotypes, this has since
% been removed

global TEMPORARYFOLDER;



if parallel==1

    timesvariant=zeros(GenomeLength,1);
    
    
    %run analysis
    parallel_params={};
    for i=1:length(SampleDirs)    
        parallel_params{end+1}={SampleDirs{i}, SampleNames{i}, ScafNames, maxFQ, onlysnps, L, TEMPORARYFOLDER};        
    end
    run_parallel_matlab_commands('generate_positions_single_sample', parallel_params, jobsubmitoptions, 1);
    
    
    %load files
    for i=1:length(SampleDirs)  
        %http://www.vsoch.com/2010/11/loading-dynamic-variables-in-a-static-workspace-in-matlab/
        pos=load([TEMPORARYFOLDER '/vcf_' SampleNames{i} '.mat']);
        if numel(pos.Positions)>2 % HC 9/13/2013
            x=chrpos2index(pos.Positions,ChrStarts);
            timesvariant(x)=timesvariant(x)+1;           
        end
        delete([TEMPORARYFOLDER '/vcf_' SampleNames{i} '.mat'])
    end

    p=find(timesvariant>0 & timesvariant <numel(SampleDirs));
    
    fprintf(['Not considering ' num2str(sum(timesvariant==numel(SampleDirs))) ' positions where all samples have a variant compared to the reference...\n'])

    
    
else

    K = 0 ; 
    Positions = zeros(L,2) ;
    Mnum = zeros(L,1) ;

    for i=1:length(SampleDirs)    
        %fprintf(1,'Strain: %g  ',i) ;    
        vcf = read_vcf_file([StrainDirs{i} '/variant.vcf']) ;     
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
                t = read_vcf_info(vcf(j).info) ;
                if t.FQ < maxFQ
                    include(K)=1;
                end
            end
        end
    end
    fprintf(1,'\n') ;  
    Positions = Positions(1:K,:) ;
    p=chrpos2index(Positions,ChrStarts);
end
