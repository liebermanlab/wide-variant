%function build_experiment_directories(SCRIPTSDIRECTORY, CLUSTERDIRECTORY, paired)

%February 2012 - November 2015
%
%Roy Kishony, Tami Lieberman, Seungsoo Kim, and Idan Yelin
% Hattie Chung: added cutadapt pre-filtering


%Takes unaligned fastq files from Illumina sequencing and produces a
%a directory for each sample containing files for downstream processing
%including alignments and variant calls


%It would be great if coverage plots were generated in this step-- they
%used to be, but the code was written very specifically for B. dolosa
%genomes and has not been generalized yet (see bottom)



%%  Important parameters that shouldn't be hardcoded in ideal version


SCRIPTSDIRECTORY='/scratch/users/mit_alm/gutevo/scripts';
REFERENCESDIRECTORY = '/scratch/users/mit_alm/gutevo/Reference_Genomes/';


paired=true;

%look for candidate indels?
find_indels=0;

%running parameters
Parallel = true ;
global RUN_ON_CLUSTER; RUN_ON_CLUSTER = 1;



%Variables that usually will not have to be changed
overwrite=0;  %overwrite
Phred_offset = 33 ; %Generally 33, was 64 in older versions. Phred_score= ascii value - Phred_offset.
adapter='CTGTCTCTTAT';

% Lung-specific variable
global LookUpScafName;
LookUpScafName = true;



%orchestra parameters
mini_q = '10';
fast_q='2:00:00';
%medium_q='6:00:00';
medium_q='12:00:00';
long_q='12:00:00';



%% Set path


path(SCRIPTSDIRECTORY,path);
fprintf(['Usings scripts directory: ' SCRIPTSDIRECTORY  '\n']);




%% Check that everything is here

SampleTable = read_sample_table ;
FilterTable = read_filter_table ;
AlignmentTable = read_alignment_table ;


%Check that alignments and filters are specified properly in FilterTable and
%Alignment Table
RefGenomes={};
for i=1:length(SampleTable)
    s=SampleTable(i) ;
    [~,ai] = ismember(s.Alignments,{AlignmentTable.Alignment}) ;
    if any(ai==0)
        error(['No alignment name found in alignment_params.csv: ' char(s(ai==0).Alignments)])
    end
    [~,fi] = ismember({AlignmentTable(ai).Filter},{FilterTable.Filter}) ;
    if any(fi==0)
        error(['No filter name found in filter_params.csv: ' char(AlignmentTable(ai(fi==0)).Filter)])
    end
    if numel(unique({FilterTable.Filter})) ~= numel({FilterTable.Filter})
        error(['Cannot have two filter settings with the same name. Check filter_params.csv'])
    end
    if numel(unique({AlignmentTable.Alignment})) ~= numel({AlignmentTable.Alignment})
        error(['Cannot have two alignment settings with the same name. Check alignment_params.csv'])
    end
    sRefGenomes={AlignmentTable(ai).Genome};
    for j=1:numel(ai)
        RefGenomes{end+1}=sRefGenomes{j};
    end
end


%Check that reference genome files are present
RefGenomes=unique(RefGenomes);
disp(RefGenomes)
for i=1:numel(RefGenomes)
    
    %check if fasta file is there and named correctly
    if ~exist([REFERENCESDIRECTORY '/' RefGenomes{i}],'dir')
        error(['Could not find a Reference Genome folder called: ' RefGenomes{i}])
    end
    if ~exist([REFERENCESDIRECTORY '/' RefGenomes{i} '/genome.fasta'],'file')
        error(['Could not find a genome.fasta file for' RefGenomes{i}])
    end
    
    %Check for presence of genebank files
    %Do not throw error-- a user may want to try many genomes
    %     before choosing a closest reference, this will save them time.
    fr = fastaread([REFERENCESDIRECTORY '/' RefGenomes{i} '/genome.fasta']) ;
    Scafs = {fr.Header} ;
    for s=1:length(Scafs)
        a=Scafs{s};
        if find(a=='|') % HC 9/16/13: proper fasta header
            f=find(a=='|',2,'last');
            fn = a(f(1)+1:f(2)-1) ;
        elseif strfind(a, 'JH') % HC 11/9/13: for smalAb55555 specific
            fn = a(1:11);
        else % HC 9/16/13: contig node
            fn = a;
        end
        
        if ~exist([REFERENCESDIRECTORY '/' RefGenomes{i} '/' fn '.gb'],'file')
            fprintf(1, ['Warning!!! could not find a genebank file called ' fn '.gb  for ' RefGenomes{i} ' \n'])
        end
    end
    
    
    %Check to see if reference genome is formatted properly for bowtie2
    if ~exist([REFERENCESDIRECTORY '/' RefGenomes{i} '/genome_bowtie2.1.bt2'],'file')
        fprintf(1,'Format reference genome for bowtie... \n') ;
        c={}; c{end+1}=['module add c3ddb/bowtie2/2.2.6; bowtie2-build -q genome.fasta genome_bowtie2'];
        d=[]; d{end+1}=[REFERENCESDIRECTORY '/' RefGenomes{i}];
        run_parallel_unix_commands_fast(c,'empty',0,d);
    end
    
end



%% Set up folders


if all([SampleTable.Lane]>0)
    %Case all lanes are not demultiplexed
    fprintf(1,'Demultiplex... \n') ;
    cmds = demultiplex(SampleTable);
    run_parallel_unix_commands_fast(cmds,medium_q,Parallel, {[pwd]});
elseif all([SampleTable.Lane]<0)
    %Case all demultiplexed already
    fprintf(1,'Copying files locally... \n') ;
    disp(SampleTable(1))
    cmds = moverenamefastqs(SampleTable);
    run_parallel_unix_commands_fast(cmds,fast_q,Parallel, {[pwd]});
else
    error('ERROR :: Currently only handles all files already demultiplexed or all files not demultiplexed. ... \n')
end

%% Trim adapter sequences

fprintf(1, 'Trim adapter sequences... \n') ; tic ;
adapter_dir = 'cutadapted';

trim_cmds = {};
trim_dirs = {};

for i = 1:length(SampleTable)
    s=SampleTable(i);
    cd(s.Sample);
      
    % check to see if next step hasn't already been processed -- if it has,
    % skip this step regardles of whether or not this folder exists.
    allprocessed=1;
    [~,ai] = ismember(s.Alignments,{AlignmentTable.Alignment});
    [~,fi] = ismember({AlignmentTable(ai).Filter},{FilterTable.Filter}) ;
    fi = unique(fi) ;
    for f=fi(:)'
        if ~exist(FilterTable(f).Filter,'dir') | ~exist([FilterTable(f).Filter '/filter_reads_1.fastq'])
            allprocessed=0;
        end
    end
    if allprocessed==1;
        cd ..
        continue; %go to next sample
    end
        
    if ~exist(adapter_dir, 'dir')    
        mkdir(adapter_dir)
    end
    
    
    % build cmds for adapter trimming
    tname_in1 = [pwd '/' s.Sample '_1.fastq'];
    tname_in2 = [pwd '/' s.Sample '_2.fastq'];
    tname_out1 = [pwd '/' adapter_dir '/cutadapt_reads_1.fastq'];
    tname_out2 = [pwd '/' adapter_dir '/cutadapt_reads_2.fastq'];
    
    if ~(exist(tname_out1,'file')) || overwrite
        trim_cmds{end+1} = ['python2.7 ~/.local/lib/python2.7/site-packages/cutadapt -a ' adapter ' ' tname_in1 ' > ' tname_out1];
        trim_dirs{end+1} = [s.Sample '/' adapter_dir];
    end
    
    if paired & (~(exist(tname_out2,'file')) || overwrite)
        trim_cmds{end+1} = ['python2.7 ~/.local/lib/python2.7/site-packages/cutadapt -a ' adapter ' ' tname_in2 ' > ' tname_out2];
        trim_dirs{end+1} = [s.Sample '/' adapter_dir];
    end
    
    cd ..
    
end

% run cutadapt with unix command
run_parallel_unix_commands_fast(trim_cmds, medium_q, Parallel, trim_dirs);
if ~exist('adapter_trimming_results','dir')
    mkdir('adapter_trimming_results')
end
prevfiles=numel(dir('adapter_trimming_results/*.sh'));
for i=1:numel(trim_cmds)
    copyfile(['run_parallel_unix_commands_fast_tmp/stdout' int2str(i) '.txt'],['adapter_trimming_results/stdout' int2str(prevfiles+i) '.txt']);
    copyfile(['run_parallel_unix_commands_fast_tmp/tmp' int2str(i) '.sh'],['adapter_trimming_results/tmp' int2str(prevfiles+i) '.sh']);
end





%% Filter reads

fprintf(1,'Filter reads... \n') ; tic ;

fparams = {};
cmds = {};
dirs ={};
for i=1:length(SampleTable)
    s=SampleTable(i) ;
    cd(s.Sample) ;
    [~,ai] = ismember(s.Alignments,{AlignmentTable.Alignment}) ;
    [~,fi] = ismember({AlignmentTable(ai).Filter},{FilterTable.Filter}) ;
    fi = unique(fi) ;
    
    for f=fi(:)'
        if ~exist(FilterTable(f).Filter,'dir')
            mkdir(FilterTable(f).Filter)
        end
        
        fname_in1=[pwd '/' adapter_dir '/cutadapt_reads_1.fastq'];
        fname_out1=[pwd '/' FilterTable(f).Filter '/filter_reads_1.fastq'];
        
        if paired %SK
            fname_in2=[pwd '/' adapter_dir '/cutadapt_reads_2.fastq'];
            fname_out2=[pwd '/' FilterTable(f).Filter '/filter_reads_2.fastq'];
            
            if ~(exist(fname_out1,'file') && exist(fname_out2,'file')) || overwrite
                % no filtering
                if strfind(FilterTable(f).Method,'nofilter') %if no filter, copy file into subdirectory-- not the best way to do this
                    fprintf(1, 'No filtering... \n');
                    cmds{end+1}=['cp ../' s.Sample '_1.fastq filter_reads_1.fastq'];
                    cmds{end+1}=['cp ../' s.Sample '_2.fastq filter_reads_2.fastq'];
                    dirs{end+1}=[s.Sample '/nofilter'];
                    dirs{end+1}=[s.Sample '/nofilter'];
                    
                    % sickle filter
                elseif strfind(FilterTable(f).Method,'sickle')
                    % modify input file to sickle
                    cmds{end+1}=['"' SCRIPTSDIRECTORY '/sickle-master/sickle" pe -f ' fname_in1 ' -r ' fname_in2 ' -t sanger -o filter_reads_1.fastq -p filter_reads_2.fastq -s singles.fastq -q ' num2str(FilterTable(f).Params(1)) ' -l ' num2str(FilterTable(f).Params(2)) '-x -n'];
                    dirs{end+1}=[s.Sample '/' FilterTable(f).Filter];
                    
                    % some other filter method
                else
                    fprintf(1, 'Some other thing in "filter"...\n');
                    fparams{end+1} =  {fname_in1, fname_in2, fname_out1, fname_out2, Phred_offset, FilterTable(f).Method, FilterTable(f).Params};
                end
            end
        else
            fprintf(1, 'Not paired! \n');
            fname_in=[pwd '/' s.Sample '.fastq'];
            fname_out=[pwd '/' FilterTable(f).Filter '/filter_reads.fastq'];
            if ~exist(fname_out1,'file') || overwrite
                if strfind(FilterTable(f).Method,'nofilter') %if no filter, copy file into subdirectory-- not the best way to do this
                    cmds{end+1}=['cp ' fname_in1 ' ' fname_out1];
                    dirs{end+1}=[s.Sample '/nofilter'];
                elseif strfind(FilterTable(f).Method,'sickle')
                    cmds{end+1}=['"' SCRIPTSDIRECTORY '/sickle-master/sickle" se -f ' fname_in1 ' -t sanger -o ' fname_out1 ' -q ' num2str(FilterTable(f).Params(1)) ' -l ' num2str(FilterTable(f).Params(2)) ' -x -n'];
                    dirs{end+1}=[s.Sample '/' FilterTable(f).Filter];
                else
                    fparams{end+1} =  {fname_in1, fname_out1, Phred_offset, FilterTable(f).Method, FilterTable(f).Params};
                end
            end
        end
    end
    cd ..
end


%run all filters that use matlab commands
if paired
    run_parallel_matlab_commands('filter_reads_paired',fparams,medium_q,Parallel);
else
    run_parallel_matlab_commands('filter_reads', fparams, medium_q, Parallel);
end

%run all filters that use unix commands
run_parallel_unix_commands_fast(cmds,medium_q,Parallel, dirs);

if ~exist('non-matlab_filter_stats','dir') & numel(cmds) > 0
    mkdir('non-matlab_filter_stats');
end

prevfiles=numel(dir('non-matlab_filter_stats/*.sh'));
for i=1:numel(cmds)
    copyfile(['run_parallel_unix_commands_fast_tmp/stdout' int2str(i) '.txt'],['non-matlab_filter_stats/stdout' int2str(prevfiles+i) '.txt']);
    copyfile(['run_parallel_unix_commands_fast_tmp/tmp' int2str(i) '.sh'],['non-matlab_filter_stats/tmp' int2str(prevfiles+i) '.sh']);
end





%% Clean up after -- if filter worked, remove cutadapted folder

for i=1:length(SampleTable)
    s=SampleTable(i) ;
    cd(s.Sample) ;
    [~,ai] = ismember(s.Alignments,{AlignmentTable.Alignment}) ;
    [~,fi] = ismember({AlignmentTable(ai).Filter},{FilterTable.Filter}) ;
    fi = unique(fi) ;
    
    for f=fi(:)'
        if exist(FilterTable(f).Filter,'dir') & exist([FilterTable(f).Filter '/filter_reads_1.fastq'])
            !rm -r cutadapted
        end
    end
    cd ..
end

%% Align


fprintf(1,'Align... \n') ;

cmds = {} ;
bowtiedirs= {} ;
all_dirs = {} ;
all_genomes = {} ;


for i=1:length(SampleTable)
    s=SampleTable(i) ;
    cd(s.Sample) ;
    for a = s.Alignments'
        
        ai = find(strcmp({AlignmentTable.Alignment},a)) ;
        ae = AlignmentTable(ai) ;
        fi = find(strcmp({FilterTable.Filter},ae.Filter)) ;
        dr = [FilterTable(fi).Filter '/' ae.Alignment] ;
        if ~exist(dr,'dir')
            mkdir(dr)
        end
        if ~exist([dr '/alignment_info'])
            save([dr '/alignment_info'], 's','ae')
        end
        
        drf=[pwd '/' ae.Filter];
        dra=[pwd '/' ae.Filter '/' ae.Alignment];
        
        if ~exist([dr '/aligned.sam']) & ~exist([dr '/aligned.sorted.bam'])
            bowtiedirs{end+1} = [REFERENCESDIRECTORY '/' ae.Genome] ;
            fprintf(['\nAligning ' s.Sample ' with ' ae.Method '!\n']);
            
            switch ae.Method
                case 'bowtie'
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie/bowtie ' ae.Param1 '--phred64-quals --max ' dra '/multialigned.fastq --un ' dra '/unaligned.fastq ' ...
                            'genome  ' drf '/filter_reads_1.fastq ' dra '/aligned.sam' ] ;
                    else
                        cmds{end+1} = ['/opt/bowtie/bowtie ' ae.Param1 ' --max ' dra '/multialigned.fastq --un ' dra '/unaligned.fastq ' ...
                            'genome  ' drf '/filter_reads_1.fastq ' dra '/aligned.sam' ] ;
                    end
                    
                case 'bowtie2'
                    if Phred_offset == 64
                        cmds{end+1} = ['bowtie2 --phred64 -x genome_bowtie2' ...
                            ' -U   ' drf '/filter_reads_1.fastq -S ' dra '/aligned.sam --un ' dra '/unaligned.fastq '];
                    else
                        cmds{end+1} = ['bowtie2 -x genome_bowtie2' ...
                            ' -U   ' drf '/filter_reads_1.fastq -S ' dra '/aligned.sam --un ' dra '/unaligned.fastq '];
                    end
                case 'bowtie2paired'
                    if Phred_offset == 64
                        cmds{end+1} = ['bowtie2  --phred64  -x genome_bowtie2' ...
                            '-1  ' drf '/filter_reads_1.fastq -2  ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam'];
                    else
                        cmds{end+1} = ['bowtie2 -x genome_bowtie2' ...
                            ' -1  ' drf '/filter_reads_1.fastq -2  ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam '];
                    end
                case 'bowtie2pairedfilter' % allows no ambiguous characters
                    if Phred_offset == 64
                        cmds{end+1} = ['bowtie2 -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --un-conc ' dra '/unaligned.fastq --phred64 -x  genome_bowtie2' ...
                            '-1  ' drf '/filter_reads_1.fastq -2  ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam'];
                    else
                        cmds{end+1} = ['bowtie2 -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --un-conc ' dra '/unaligned.fastq -x  genome_bowtie2' ...
                            ' -1 ' drf '/filter_reads_1.fastq -2  ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam '];
                    end
                case 'bowtie2pairedxt' % allows no ambiguous characters
                    if Phred_offset == 64
                        cmds{end+1} = ['bowtie2 -X 2000 --no-mixed --dovetail --very-sensitive --n-ceil 0,0.01 --un-conc ' dra '/unaligned.fastq --phred64 -x  genome_bowtie2' ...
                            '-1 ' drf '/filter_reads_1.fastq -2 ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam'];
                    else
                        cmds{end+1} = ['bowtie2 -X 2000 --no-mixed --dovetail --very-sensitive --n-ceil 0,0.01 --un-conc ' dra '/unaligned.fastq -x  genome_bowtie2' ...
                            ' -1 ' drf  '/filter_reads_1.fastq -2 ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam '];
                    end
                    
                    
            end
        end
        all_dirs{end+1} = [pwd '/' dr] ;
        all_genomes{end+1} = ae.Genome ;
    end
    cd ..
end

run_parallel_unix_commands_fast(cmds, long_q, Parallel, bowtiedirs);

%save results
if ~exist('alignment_stats','dir') & numel(cmds) > 0
    !mv run_parallel_unix_commands_fast_tmp alignment_stats
end




 

%this will only work if bowtie2 was used. creates a data strucutre for
%later visualization of alignments
%summarize_alignments



%% compress and sort

fprintf(1,'Compress to .bam and sort ... \n') ;

dirs = {} ;
cmds = {} ;
for i=1:length(all_dirs)
    if ~exist([all_dirs{i} '/aligned.sorted.bam'],'file') ;
        dirs{end+1} = all_dirs{i} ;
    end
end

run_parallel_unix_commands_fast({'samtools view -bS -o aligned.bam aligned.sam; samtools sort aligned.bam aligned.sorted; rm aligned.sam; rm aligned.bam'},long_q,Parallel,dirs);





%% Make pileup files for position-specific statistics

fprintf(1,'Pileup... \n') ;

dirs = {} ;
cmds = {} ;
for i=1:length(all_dirs)
    if ~exist([all_dirs{i} '/strain.pileup'],'file') ; %remove pileup2
        %option -s -O increases information output (mapping quality and
        %position on read
        %-B disables BAQ computation
        %-d sets maximum read depth
        cmds{end+1} = ['samtools mpileup -q30 -x -s -O -d3000 -f "' REFERENCESDIRECTORY '/' all_genomes{i} '/genome.fasta" aligned.sorted.bam > strain.pileup'] ;
        dirs{end+1} = all_dirs{i} ;
    end
end

run_parallel_unix_commands_fast(cmds,medium_q,Parallel,dirs);



%% Make vcf files for consensus calling

dirs = {} ;
cmds = {} ;
for i=1:length(all_dirs)
    if ~exist([all_dirs{i} '/variant.vcf'],'file') ;
        %HC 11/6/2013: -B disables BAQ computation, consistent with pileup
        %generation   
       cmds{end+1} = ['samtools mpileup -q30 -t SP -d3000 -uvf"' REFERENCESDIRECTORY '/' all_genomes{i} '/genome.fasta" aligned.sorted.bam > strain; bcftools call -c -Ov -o strain.vcf strain; bcftools view -v snps -q .75 strain.vcf > variant.vcf'] ;
       dirs{end+1} = all_dirs{i} ;
    end
end

run_parallel_unix_commands_fast(cmds,medium_q,Parallel,dirs);


%%  Store FQ scores from strain.vcf into quals.mat


for i=1:length(all_dirs)
    if ~exist([all_dirs{i} '/quals.mat'],'file') ;
        cmds{end+1} = ['matlab -r "path(' char(39) SCRIPTSDIRECTORY char(39) ',path); vcf_to_quals(' char(39) all_dirs{i} char(39) ',' char(39) REFERENCESDIRECTORY '/' all_genomes{i} char(39) ')"'];
    end
end

run_parallel_unix_commands_fast(cmds,medium_q,Parallel,{'.'});


%% Call indels using dindel

if find_indels==1
    
    fprintf(1,'Finding candidate indels...\n');
    
    dirs = {} ;
    cmds={};
    for i=1:length(all_dirs)
        if ~exist([all_dirs{i} '/sample.dindel_output.variants.txt'],'file') ;
            cmds{end+1} = ['dindel --analysis getCIGARindels --bamFile aligned.sorted.bam --outputFile sample.dindel_output' ...
                ' --ref ' ref_folder '/Reference_Genomes/' all_genomes{i} '/genome.fasta'];
            dirs{end+1} = all_dirs{i} ;
        end
    end
    
    run_parallel_unix_commands_fast(cmds,long_q,Parallel,dirs);
    
end

%% Create diversity.mat (from pileup files) for diversity calling

fprintf(1,'Create diversity.mat... \n') ;

dirs = {} ;
for i=1:length(all_dirs)
    if ~exist([all_dirs{i} '/diversity.mat'],'file') ;
        dirs{end+1} = all_dirs{i} ;    
    end
end

run_parallel_unix_commands_fast({['matlab -r "path(' char(39) SCRIPTSDIRECTORY char(39) ',path); pileup_to_diversity_matrix"']},long_q,Parallel,dirs);



fprintf(1,'\nBuilding .bai files for viewing alignment... \n') ;

%% Rename bam first, then make bai files, then delete temporary files


cmds_bai = {} ;
dirs = {} ;

for i=1:length(SampleTable)
    s=SampleTable(i) ;
    for a = s.Alignments'
        ai = find(strcmp({AlignmentTable.Alignment},a)) ;
        ae = AlignmentTable(ai) ;
        fi = find(strcmp({FilterTable.Filter},ae.Filter)) ;
        dr = [s.Sample '/' FilterTable(fi).Filter '/' ae.Alignment] ;
        if ~exist([dr '/aligned.sorted.bam.bai'])
            dirs{end+1} = [pwd '/' dr] ;
        end
    end
end

run_parallel_unix_commands_fast({'samtools index aligned.sorted.bam'},long_q,Parallel,dirs);

%% Clean up after yourself

for i=1:length(all_dirs)
    cd(all_dirs{i})
    d=dir('diversity.mat');
    q=dir('quals.mat');
    if ~isempty(d) & ~isempty(q) & d.bytes > 0 & q.bytes > 0
        !rm strain
        !rm strain.pileup
        !rm strain.vcf
    else
        error(['Error: quals.mat and diversity mat not generated successfully for directory ' all_dirs{i}]);
    end
end
