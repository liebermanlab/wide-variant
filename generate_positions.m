function p = generate_positions(SampleDirs, SampleNames, ScafNames, ChrStarts, L, parallel, jobsubmitoptions)      
% unlike generate_mutations_genotypes, this functions creates a list of mutated positions (chromosome, position) rather than a list of mutaions (chromosome, position, ref, alt).
% changes are mainly in not reading reference and alelle columns.  

% an earlier version created a structure called Genotypes, this has since
% been removed




if parallel==1

    
    
    
    %run analysis
    parallel_params={};
    for i=1:length(SampleDirs)    
        parallel_params{end+1}={SampleDirs{i}, SampleNames{i}, ScafNames, L};        
    end
    run_parallel_matlab_commands('generate_positions_single_sample', parallel_params, jobsubmitoptions, 1);
    
    
    %load files
    p=[];
    for i=1:length(SampleDirs)    
        %http://www.vsoch.com/2010/11/loading-dynamic-variables-in-a-static-workspace-in-matlab/
        pos=load(['vcf_' SampleNames{i} '.mat']);
        p=union(p, chrpos2index(pos.Positions,ChrStarts));
        delete(['vcf_' SampleNames{i} '.mat'])
    end

   
        
    
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
            end
        end
    end
    fprintf(1,'\n') ;  
    Positions = Positions(1:K,:) ;
    p=chrpos2index(Positions,ChrStarts);
end
