function run_parallel_matlab_commands(func_name, params, options, Parallel)

global SCRIPTSPATH

if isempty(params)
    return
end


func = str2func(func_name) ;
if Parallel
    jm=findResource('scheduler','type','lsf');
    set(jm, 'ClusterMatlabRoot', '/opt/matlab')
    job=createJob(jm, 'PathDependencies', {SCRIPTSPATH}); %'FileDependencies', {'salkja'}
    set(jm, 'SubmitArguments',['-R "rusage[matlab_dc_lic=1]" -q ' options] );
    
    %fprintf('Preparing parallel matlab jobs...\n')
    for i=1:length(params)
        createTask(job, func, 0, params{i});
    end
    
    
    
    submit(job)
    fprintf(['Sent ' num2str(length(params)) ' jobs. Waiting...\n'])
    waitForState(job)
    
    for i=1:length(params)
        if (isempty(job.Tasks(i).FinishTime)) | (~isempty(job.Tasks(i).ErrorMessage))
            if isempty(job.Tasks(i).FinishTime)
                error(['Matlab job [' num2str(i) '] failed because of time limitation'])
            else
                disp(job.Task(i).ErrorMessage)
                error(['Matlab job [' num2str(i) '] failed (others may have failed also)'])
            end
        end
    end
    
    
else
    for i=1:length(params)
        func(params{i}{:});
    end 
end
    
end


