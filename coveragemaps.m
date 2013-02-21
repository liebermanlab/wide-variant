% coveragemaps.m
%
% run coverage_map on all Staph samples
% 
% Seungsoo Kim
% December 30, 2012


IsolateTable = read_isolate_table ;
Filter = 'trim_end25_20';
Alignment = 'bowtie2';

params = {};
for i = 1:length(IsolateTable)
    s=IsolateTable(i);
    if ~exist([s.Sample '/' Filter '/' Alignment '/ coverage_map.tif'],'file')
        params{end+1} = {[s.Sample '/' Filter '/' Alignment]};
    end
end

run_parallel_matlab_commands('coverage_map',params,'sysbio_12h',0);