function run_parallel_unix_commands(cmds, qname) 

jm=findResource('scheduler','type','lsf');
set(jm, 'ClusterMatlabRoot', '/opt/matlab')
job=createJob(jm, 'PathDependencies', {'/files/SysBio/KISHONY LAB/Illumina_pipeline/scripts'}); %'FileDependencies', {'salkja'}
 set(jm, 'SubmitArguments',['-R "rusage[matlab_dc_lic=' num2str(length(cmds)) ']" -q ' qname] );

for i=1:length(cmds)
    createTask(job, @eval, 0, {['! ' cmds{i}]});
end    

submit(job)
waitForState(job)
destroy(job) %supresses writing of error files
end

