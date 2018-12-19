function send_jobs_to_cluster(cmds, maxtime, Parallel, dirs)

%Edited November 2015 to work on c3ddb
%Edited 2018.10.31 (Arolyn): Lines 40 & 92 to dispatch jobs to defq or
%sched_mem1TB, whichever is available first

% function send_jobs_to_cluster(cmds, qname)
% a function that runs unix commands (given in a cell array cmds) in parallel by generating temp sh files
% that also produce output files and waiitng for all output files to
% appear. uses queue given in qname. Can run on mulitple dirs in a ingle
% command, Or multiple commands signle directory.

% edited by Seungsoo Kim, December 17, 2012
% to search for output line including 'Subject: ' and then check whether
% 'Done'

%Comes in 3 flavors, designated by Parallel
%Parallel = 1 -> run on cluster, wait for output
%Parallel = 2 -> run on cluster, do not wait for output
%Otherwise-- run sequentially, wait for output



if isempty(cmds) | isempty(dirs)
    return
end



if Parallel==1
    dr = 'temp_batch_files_and_output' ;
    outs={}; %list of orchestra output files to expect.
    if ~exist(dr,'dir')
        mkdir(dr) ;
    end
    delete([dr ,'/*']) ;
    for i=1:max(length(cmds),length(dirs)) 
        fname = [dr '/tmp' num2str(i) '.sh'] ;
        fid = fopen(fname,'w') ;
        fprintf(fid,'#!/bin/bash\n') ;
        fprintf(fid,'#SBATCH -p defq,sched_mem1TB\n') ;
        fprintf(fid,'#SBATCH -n 1\n') ;
        fprintf(fid,'#SBATCH --time=%s\n',maxtime) ;
        fprintf(fid,'#SBATCH --mem=50000\n') ;
        outs{end+1}=[dr '/stdout' num2str(i) '.txt'];
        fprintf(fid,'#SBATCH -o %s\n',outs{end}) ;
        fprintf(fid,'#SBATCH -e %s/stderr%g.txt\n',dr,i) ;
        fprintf(fid,'module add c3ddb/bamtools/2.4.0\n') ;
        fprintf(fid,'module add c3ddb/bcftools/1.2\n') ;
        fprintf(fid,'module add c3ddb/bowtie2/2.2.6 \n') ;
        fprintf(fid,'module add c3ddb/htslib/1.2.1\n') ;
        fprintf(fid,'module add c3ddb/samtools/1.2\n') ;
        fprintf(fid,'module add c3ddb/sickle/1.33\n') ;
        fprintf(fid,'module add c3ddb/samtools/1.2\n') ;
        fprintf(fid,'module add c3ddb/python/2.7.11\n') ;
        fprintf(fid,'pip install --user --upgrade cutadapt\n') ;
        fprintf(fid,'pip install --user FileUtils\n') ;
        fprintf(fid,'module add mit/matlab/2015b\n') ;
        fprintf(fid,'module add c3ddb/deML/2016-05-25\n') ;
        fprintf(fid,'cd "%s"\n',dirs{min(i,end)}) ;
        fprintf(fid,'%s\n',cmds{min(i,end)}) ;
        fprintf(fid,'echo Done!!!') ;
        fclose(fid) ;
        eval(sprintf('!sbatch %s',fname))
    end
    
    done = 0 ;
    numjobs=max(length(cmds),length(dirs));
    while done < numjobs
        pause(30) ;
        for i=1:length(outs)
            if exist(outs{i},'file')
                fid=fopen(outs{i});
                l=fgetl(fid);
                while sum(l==-1)==0
                    if strfind(l,'Done!!!')
                        done=done+1;
                        outs{i}=[];
                    end
                    l=fgetl(fid);
                end
                fclose(fid);
            end    
        end
           fprintf(1,['Waiting on ' num2str(numjobs - done) 'jobs \n']);
    end
    
elseif Parallel==2
    for i=1:max(length(cmds),length(dirs))
        fname = ['tmp' num2str(i) '.sh'] ;
        fid = fopen(fname,'w') ;
        fprintf(fid,'#!/bin/bashn') ;
        fprintf(fid,'#SBATCH -p defq,sched_mem1TB\n') ;
        fprintf(fid,'#SBATCH -n 1\n') ;
        fprintf(fid,'#SBATCH --time=%s\n',maxtime) ; 
        fprintf(fid,'#SBATCH -o stdout%g.txt\n',i) ; 
        fprintf(fid,'#SBATCH -e stderr%g.txt\n',i) ; 
        fprintf(fid,'module add c3ddb/bamtools/2.4.0\n') ;
        fprintf(fid,'module add c3ddb/bcftools/1.2\n') ;
        fprintf(fid,'module add c3ddb/bowtie2/2.2.6 \n') ;
        fprintf(fid,'module add c3ddb/htslib/1.2.1\n') ;
        fprintf(fid,'module add c3ddb/samtools/1.3\n') ;
        fprintf(fid,'module add c3ddb/sickle/1.33\n') ;
        fprintf(fid,'module add c3ddb/samtools/1.2\n') ;
        fprintf(fid,'module add c3ddb/python/2.7.10\n') ;
        fprintf(fid,'pip install --user cutadapt==1.9.1\n') ;
        fprintf(fid,'module add mit/matlab/2015b\n') ;        
        fprintf(fid,'cd "%s"\n',dirs{min(i,end)}) ; 
        fprintf(fid,'%s\n',cmds{min(i,end)}) ;
        fclose(fid) ;
        eval(sprintf('!sbatch %s',fname))
    end
    fprintf('Continuing without waiting for last batch of jobs to finish...\n')
else
    cdr = pwd ;
    for i=1:max(length(cmds),length(dirs))
        cd(dirs{min(i,end)})
        eval(['!' cmds{min(i,end)}]) ;
        cd(cdr)
    end
end


end

