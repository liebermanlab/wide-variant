function build_experiment_directories(scriptspath)

%February 2012 - March 2013
%
%Roy Kishony, Tami Lieberman, Seungsoo Kim, and Idan Yelin


%Takes unaligned fastq files from Illumina sequencing and produces a
%a directory for each sample containing files for downstream processing
%including alignments and variant calls


%It would be great if coverage plots were generated in this step-- they
%used to be, but the code was written very specifically for B. dolosa
%genomes and has not been generalized yet (see bottom)



%%  Important parameters that shouldn't be hardcoded in ideal version


%running parameters
paired = true; % could be acquired from files
Parallel = true ;
global RUN_ON_CLUSTER; RUN_ON_CLUSTER = 1;

%orchestra parameters
fast_q='sysbio_15m';
demultiplex_q='sysbio_2h';
filter_q='sysbio_2h';
alignment_q='sysbio_12h';
processing_q='sysbio_12h';
plotting_q='sysbio_int_2h';

%Variables that usually will not have to be changed
overwrite=0;  %overwrite
Phred_offset = 33 ; %Generally 33, was 64 in older versions. Phred_score= ascii value - Phred_offset.


%% Set path


global SCRIPTSPATH;


if nargin > 0
    SCRIPTSPATH = scriptspath;
else
   a=pwd; SCRIPTSPATH=[a(1:find(double(a)==47,1,'last')) 'scripts']; 
   if ~exist(SCRIPTSPATH, 'dir')
       fprintf(['Could not find ' SCRIPTSPATH '...\n'])
       error('Error: Must run from an experiment folder, with scripts folder in parent directory')
   else
       path(SCRIPTSPATH,path);
   end
end

fprintf(['Usings scripts directory: ' SCRIPTSPATH  '\n']);



%% Set main folder
if RUN_ON_CLUSTER == 1
    mainfolder='/files/SysBio/KISHONY LAB/illumina_pipeline';
else
    mainfolder='/Volumes/sysbio/KISHONY LAB/illumina_pipeline';
end


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
for i=1:numel(RefGenomes)
    
    %check if fasta file is there and named correctly
    if ~exist([mainfolder '/Reference_Genomes/' RefGenomes{i}],'dir')
        error(['Could not find a Reference Genome folder called: ' RefGenomes{i}])
    end
    if ~exist([mainfolder '/Reference_Genomes/' RefGenomes{i} '/genome.fasta'],'file')
        error(['Could not find a genome.fasta file for' RefGenomes{i}])
    end
    
    %Check for presence of genebank files
    %Do not throw error-- a user may want to try many genomes
    %     before choosing a closest reference, this will save them time.
    fr = fastaread([mainfolder '/Reference_Genomes/' RefGenomes{i} '/genome.fasta']) ;
    Scafs = {fr.Header} ;
    for s=1:length(Scafs)
        a=Scafs{s};
        f=find(a=='|',2,'last');fn = a(f(1)+1:f(2)-1) ;
        if ~exist([mainfolder '/Reference_Genomes/' RefGenomes{i} '/' fn '.gb'],'file')
            fprintf(1, ['Warning!!! could not find a genebank file called ' fn '.gb  for ' RefGenomes{i} ' \n'])
        end
    end
    
    
    %Check to see if reference genome is formatted properly for bowtie2
    if ~exist([mainfolder '/Reference_Genomes/' RefGenomes{i} '/genome_bowtie2.1.bt2'],'file')
        fprintf(1,'Format reference genome for bowtie... \n') ; tic ;
        c={}; c{end+1}=['/opt/bowtie2/bowtie2-build -q genome.fasta genome_bowtie2'];
        d=[]; d{end+1}=[mainfolder '/Reference_Genomes/' RefGenomes{i}];
        run_parallel_unix_commands_fast(c,'empty',0,d);
    end
    


end





%% Set up folders


if all([SampleTable.Lane]>0)
    %Case all lanes are not demultiplexed
    fprintf(1,'Demultiplex... \n') ; tic ;
    cmds = demultiplex(SampleTable);
    run_parallel_unix_commands_fast(cmds,demultiplex_q,Parallel, {[pwd]});
elseif all([SampleTable.Lane]<0)
    %Case all demultiplexed already
    fprintf(1,'Copying files locally... \n') ; tic ;
    cmds = moverenamefastqs(SampleTable);
    run_parallel_unix_commands_fast(cmds,demultiplex_q,Parallel, {[pwd]});
else
    error('ERROR :: Currently only handles all files already demultiplexed or all files not demultiplexed. ... \n')
end

%% Trim adapter sequences 

fprintf(1, 'Trim adapter sequences... \n') ; tic ; 

tparams = {}; 
cmds = {}; 
dirs = {}; 

for i = 1:length(SampleTable)
    s=SampleTable(i);
    cd(s.Sample);
    
end

%% Filter reads

fprintf(1,'Filter reads... \n') ; tic ;

fparams = {};
cmds = {} ;
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
               
        
        if paired %SK
            fname_in1=[pwd '/' s.Sample '_1.fastq'];
            fname_in2=[pwd '/' s.Sample '_2.fastq'];
            fname_out1=[pwd '/' FilterTable(f).Filter '/filter_reads_1.fastq'];
            fname_out2=[pwd '/' FilterTable(f).Filter '/filter_reads_2.fastq'];
            if ~(exist(fname_out1,'file') && exist(fname_out2,'file')) || overwrite
                if strfind(FilterTable(f).Method,'nofilter') %if no filter, copy file into subdirectory-- not the best way to do this
                    cmds{end+1}=['cp ../' s.Sample '_1.fastq filter_reads_1.fastq'];
                    cmds{end+1}=['cp ../' s.Sample '_2.fastq filter_reads_2.fastq'];
                    dirs{end+1}=[s.Sample '/nofilter'];
                    dirs{end+1}=[s.Sample '/nofilter'];
                elseif strfind(FilterTable(f).Method,'sickle')
                    cmds{end+1}=['"' SCRIPTSPATH '/sickle-master/sickle" pe -f ../' s.Sample '_1.fastq -r ../' s.Sample '_2.fastq -t sanger -o filter_reads_1.fastq -p filter_reads_2.fastq -s singles.fastq -q ' num2str(FilterTable(f).Params(1)) ' -l ' num2str(FilterTable(f).Params(2)) '-x -n'];
                    dirs{end+1}=[s.Sample '/' FilterTable(f).Filter];
                else
                    fparams{end+1} =  {fname_in1, fname_in2, fname_out1, fname_out2, Phred_offset, FilterTable(f).Method, FilterTable(f).Params};
                end
            end
        else
            fname_in=[pwd '/' s.Sample '.fastq'];
            fname_out=[pwd '/' FilterTable(f).Filter '/filter_reads.fastq'];
            if ~exist(fname_out,'file') || overwrite
                if strfind(FilterTable(f).Method,'nofilter') %if no filter, copy file into subdirectory-- not the best way to do this
                    cmds{end+1}=['cp ../' s.Sample '.fastq filter_reads.fastq'];
                    dirs{end+1}=[s.Sample '/nofilter'];
                elseif strfind(FilterTable(f).Method,'sickle')
                    cmds{end+1}=['"' SCRIPTSPATH '/sickle-master/sickle" se -f ../' s.Sample '.fastq -t sanger -o filter_reads.fastq -q ' num2str(FilterTable(f).Params(1)) ' -l ' num2str(FilterTable(f).Params(2)) ' -x -n'];
                    dirs{end+1}=[s.Sample '/' FilterTable(f).Filter];
                else
                    fparams{end+1} =  {[pwd '/' s.Sample '.fastq'], [pwd '/' FilterTable(f).Filter '/filter_reads.fastq'], Phred_offset, FilterTable(f).Method, FilterTable(f).Params};
                end
            end
        end
    end
    cd ..
end


%run all filters that use matlab commands
if paired
    run_parallel_matlab_commands('filter_reads_paired',fparams,filter_q,Parallel);
else
    run_parallel_matlab_commands('filter_reads', fparams, filter_q, Parallel);
end

%run all filters that use unix commands
run_parallel_unix_commands_fast(cmds,filter_q,Parallel, dirs);

if ~exist('non-matlab_filter_stats','dir') & numel(cmds) > 0
    mkdir('non-matlab_filter_stats')
    !mv run_parallel_unix_commands_fast_tmp/*  non-matlab_filter_stats/
end



%% Align


fprintf(1,'Align... \n') ; tic ;

cmds = {} ;
bowtiedirs= {} ; 
dirs = {} ;
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
        
        if ~exist([dr '/aligned.sorted.bam'])
            
            dirs{end+1} = [pwd '/' dr] ;
            bowtiedirs{end+1} = [mainfolder '/Reference_Genomes/' ae.Genome] ;
            
            switch ae.Method
                case 'bowtie'
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie/bowtie ' ae.Param1 '--phred64-quals --max ' dra '/multialigned.fastq --un ' dra '/unaligned.fastq ' ...
                            'genome  ' drf '/filter_reads.fastq ' dra '/aligned.sam' ] ;
                    else
                        cmds{end+1} = ['/opt/bowtie/bowtie ' ae.Param1 ' --max ' dra '/multialigned.fastq --un ' dra '/unaligned.fastq ' ...
                            'genome  ' drf '/filter_reads.fastq ' dra '/aligned.sam' ] ;
                    end
                    
                case 'bowtie2'
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 --phred64 -x genome_bowtie2' ...
                            ' -U   ' drf '/filter_reads.fastq -S ' dra '/aligned.sam --un ' dra '/unaligned.fastq '];
                    else
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -x genome_bowtie2' ...
                            ' -U   ' drf '/filter_reads.fastq -S ' dra '/aligned.sam --un ' dra '/unaligned.fastq '];
                    end
                case 'bowtie2paired'
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie2/bowtie2  --phred64  -x genome_bowtie2' ...
                            '-1  ' drf '/filter_reads_1.fastq -2  ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam'];
                    else
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -x genome_bowtie2' ...
                            ' -1  ' drf '/filter_reads_1.fastq -2  ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam '];
                    end
                case 'bowtie2pairedfilter' % allows no ambiguous characters
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --un-conc ' dra '/unaligned.fastq --phred64 -x  genome_bowtie2' ...
                            '-1  ' drf '/filter_reads_1.fastq -2  ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam'];
                    else
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --un-conc ' dra '/unaligned.fastq -x  genome_bowtie2' ...
                            ' -1 ' drf '/filter_reads_1.fastq -2  ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam '];
                    end
                case 'bowtie2pairedxt' % allows no ambiguous characters
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -X 2000 --no-mixed --dovetail --very-sensitive --n-ceil 0,0.01 --un-conc ' dra '/unaligned.fastq --phred64 -x  genome_bowtie2' ...
                            '-1 ' drf '/filter_reads_1.fastq -2 ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam'];
                    else
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -X 2000 --no-mixed --dovetail --very-sensitive --n-ceil 0,0.01 --un-conc ' dra '/unaligned.fastq -x  genome_bowtie2' ...
                            ' -1 ' drf  '/filter_reads_1.fastq -2 ' drf '/filter_reads_2.fastq -S ' dra '/aligned.sam '];
                    end
                    
                    
            end
        end
        all_dirs{end+1} = [pwd '/' dr] ;
        all_genomes{end+1} = ae.Genome ;
    end
    cd ..
end


if ~exist('alignment_stats','dir')
    run_parallel_unix_commands_fast(cmds,alignment_q,Parallel,bowtiedirs);
    mkdir('alignment_stats')
    !mv run_parallel_unix_commands_fast_tmp/*  alignment_stats/
end

%this will only work if bowtie two was used. creates a data strucutre for
%later visualization of alignments
%summarize_alignments 

fprintf(1,'Create .bam and sort... \n') ;
run_parallel_unix_commands_fast({'/opt/samtools/bin/samtools view -bS -o aligned.bam aligned.sam'},alignment_q,Parallel,dirs);
run_parallel_unix_commands_fast({'/opt/samtools/bin/samtools sort aligned.bam aligned.sorted'},alignment_q,Parallel,dirs);
run_parallel_unix_commands_fast({'rm aligned.sam'},alignment_q,0,dirs);
run_parallel_unix_commands_fast({'rm aligned.bam'},alignment_q,0,dirs);


%% Make pileup files for position-specific statistics 

fprintf(1,'Pileup... \n') ; tic ;

dirs = {} ;
cmds = {} ;
for i=1:length(all_dirs)
    if ~exist([all_dirs{i} '/strain.pileup'],'file') ; %remove pileup2
        %option -s -O increases information output (mapping quality and
        %position on read
        %-B disables BAQ computation
        %-d sets maximum read depth
        if Phred_offset==64
            cmds{end+1} = ['/opt/samtools/bin/samtools mpileup -q30 -6 -O -s -B -d3000 -f  "' mainfolder '/Reference_Genomes/' all_genomes{i} '/genome.fasta" aligned.sorted.bam > strain.pileup'] ;
        else
            cmds{end+1} = ['/opt/samtools/bin/samtools mpileup -q30 -s -O -B -d3000 -f "' mainfolder '/Reference_Genomes/' all_genomes{i} '/genome.fasta" aligned.sorted.bam > strain.pileup'] ;
        end
        dirs{end+1} = all_dirs{i} ;
    end
end

run_parallel_unix_commands_fast(cmds,processing_q,Parallel,dirs);



%% Make vcf files for consensus calling



dirs = {} ;
cmds = {} ;
for i=1:length(all_dirs)
    if ~exist([all_dirs{i} '/variant.vcf'],'file') ;
        if Phred_offset==64
            cmds{end+1} = ['/opt/samtools/bin/samtools mpileup -q30 -6 -S -ugf "' mainfolder '/Reference_Genomes/' all_genomes{i} '/genome.fasta" aligned.sorted.bam > strain'] ;
        else
            cmds{end+1} = ['/opt/samtools/bin/samtools mpileup -q30 -S -ugf "' mainfolder '/Reference_Genomes/' all_genomes{i} '/genome.fasta" aligned.sorted.bam > strain'] ;
        end
        
        dirs{end+1} = all_dirs{i} ;
    end
end

run_parallel_unix_commands_fast(cmds,processing_q,Parallel,dirs);


%vcf files are called for the purpose of making calls of difference from
%reference
run_parallel_unix_commands_fast({'/opt/samtools/bin/bcftools view -g strain > strain.vcf'},alignment_q,Parallel,dirs);
run_parallel_unix_commands_fast({'/opt/samtools/bin/bcftools view -vS strain.vcf > variant.vcf'},alignment_q,Parallel,dirs);




%% Create diversity.mat (from pileup files) for diversity calling

fprintf(1,'Create diversity.mat... \n') ; tic ;


params = {} ;
for i=1:length(all_dirs)
    if ~exist([all_dirs{i} '/diversity.mat'],'file') ;
        params{end+1} =  {[all_dirs{i} '/strain.pileup'], [all_dirs{i} '/diversity.mat'], all_genomes{i}};
    end
end

run_parallel_matlab_commands('pileup_to_diversity_matrix', params, processing_q, 1);



fprintf(1,'Building .bai files for viewing alignment... \n') ; tic ;


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

run_parallel_unix_commands_fast({'/opt/samtools/bin/samtools index aligned.sorted.bam'},alignment_q,Parallel,dirs);



end


