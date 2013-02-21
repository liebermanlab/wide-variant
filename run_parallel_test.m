function run_parallel_test(qname) 

filenames={
    'test1.txt'
    'test2.txt'
    'test3.txt'
 };


jm=findResource('scheduler','type','lsf');
set(jm, 'ClusterMatlabRoot', '/opt/matlab')
job=createJob(jm, 'PathDependencies', {'/files/SysBio/KISHONY LAB/Illumina_pipeline/scripts'}); %'FileDependencies', {'salkja'}
 set(jm, 'SubmitArguments',['-R "rusage[matlab_dc_lic=' num2str(length(filenames)) ']" -q ' qname] );


for i=1:length(filenames)
    createTask(job, @test, 0, {filenames{i}});
end    

submit(job)
waitForState(job)
destroy(job) %supresses writing of error files
end

