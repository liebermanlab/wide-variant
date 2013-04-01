function run_parallel_unix_commands_fast(cmds, qname, Parallel, dirs)
% function run_parallel_unix_commands_fast(cmds, qname)
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

if nargin < 4
    dirs = {'.'} ;
end

if isempty(cmds)
    return
end

if isempty(dirs)
    return
end





if Parallel==1
    dr = 'run_parallel_unix_commands_fast_tmp' ;
    outs={}; %list of orchestra output files to expect.
    if ~exist(dr,'dir')
        mkdir(dr) ;
    end
    delete([dr ,'/*']) ;
    for i=1:max(length(cmds),length(dirs))
        fname = sprintf('%s/tmp%g.sh',dr,i) ;
        oname = sprintf('%s/out%g.txt',dr,i) ;
        outs{end+1}=oname;
        fid = fopen(fname,'w') ;
        fprintf(fid,'cd "%s"\n',dirs{min(i,end)}) ;
        fprintf(fid,'%s\n',cmds{min(i,end)}) ;
        fclose(fid) ;
        eval(sprintf('!chmod +x %s',fname))
       % eval(sprintf('!bsub -q %s ./%s',qname,fname))
        eval(sprintf('!bsub -q %s -o %s ./%s',qname,oname,fname))
    end
    
    done = 0 ;
    while done < max(length(cmds),length(dirs))
        pause(10) ;
        for i=1:length(outs)
            if exist(outs{i})
                fid=fopen(outs{i});
                l=fgetl(fid);
                while isempty(strfind(l,'Subject: Job '))
                    l=fgetl(fid);
                end
                if strfind(l,'Done')
                    done=done+1;
                    outs{i}=[];
                else
                    error('A job failed. Check run_parallel_unix_commands_fast_tmp for error message')
                end
            end    
        end
        disp(done) ;
    end
    
elseif Parallel==2
    for i=1:max(length(cmds),length(dirs))
        fname = ['tmp' num2str(i) '.sh'] ;
        oname = ['out' num2str(i) '.txt'] ;
        fid = fopen(fname,'w') ;
        fprintf(fid,'cd "%s"\n',dirs{min(i,end)}) ;
        fprintf(fid,'%s\n',cmds{min(i,end)}) ;
        fclose(fid) ;
        eval(sprintf('!chmod +x %s',fname))
        eval(sprintf('!bsub -q %s -o %s ./%s',qname,oname,fname))
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

