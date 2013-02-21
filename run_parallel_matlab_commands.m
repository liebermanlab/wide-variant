function run_parallel_matlab_commands(func_name, params, options, Parallel)


if isempty(params)
    return
end


func = str2func(func_name) ;
if Parallel
    jm=findResource('scheduler','type','lsf');
    set(jm, 'ClusterMatlabRoot', '/opt/matlab')
    job=createJob(jm, 'PathDependencies', {'/files/SysBio/KISHONY LAB/illumina_pipeline/scripts'}); %'FileDependencies', {'salkja'}
    set(jm, 'SubmitArguments',['-R "rusage[matlab_dc_lic=1]" -q ' options] );
    
    %fprintf('Preparing parallel matlab jobs...\n')
    for i=1:length(params)
        createTask(job, func, 0, params{i});
    end
    
    submit(job)
    fprintf(['Sent ' num2str(length(params)) ' jobs. Waiting...\n'])
    waitForState(job)
    %destroy(job) %supresses writing of error files
else
    for i=1:length(params)
        func(params{i}{:});
    end 
end
    
end


