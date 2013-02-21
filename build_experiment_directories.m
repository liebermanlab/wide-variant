
%Feb 2012 - Feb 2013
%
%Roy Kishony, Tami Lieberman, Seungsoo Kim, and Idan Yelin
%Takes unaligned fastq files from Illumina sequencing and produces a
%a directory for each sample containing files for downstream processing
%including alignments and variant calls


%Make sure to that dependencies in run_parallel_matlab_commands are
%properly specified

%Either files are all already demultiplexed, or none are

%It would be great if coverage plots were generated in this step-- they
%used to be, but the code was written very specifically for B. dolosa
%genomes and has not been generalized yet (see bottom)

%-------------------------------------------------------------------------
%   Important parameters that shouldn't be hardcoded in ideal version


%This could/should be taken from fastq files, 64 in older versions
Phred_offset = 33 ; 

paired = true; 
Parallel = true ;

fast_q='sysbio_15m';
demultiplex_q='sysbio_2h';
filter_q='sysbio_2h';
alignment_q='sysbio_12h';
processing_q='sysbio_2h';
plotting_q='sysbio_int_2h';

%overwrite
overwrite=0;



%-------------------------------------------------------------------------

global RUN_ON_CLUSTER
RUN_ON_CLUSTER = 1;

path('/files/SysBio/KISHONY LAB/illumina_pipeline/scripts/',path);



%-------------------------------------------------------------------------
% Set up folders

SampleTable = read_sample_table ;
FilterTable = read_filter_table ;
AlignmentTable = read_alignment_table ;


if all([SampleTable.Lane]>0)
    %Case all lanes are not demultiplexed
    fprintf(1,'Demultiplex... \n') ; tic ; 
    cmds = demultiplex(SampleTable);
    run_parallel_unix_commands_fast(cmds,demultiplex_q,Parallel);
elseif all([SampleTable.Lane]<0)
    %Case all demultiplexed already
    fprintf(1,'Copying files locally... \n') ; tic ;
    cmds = moverenamefastqs(SampleTable);
    run_parallel_unix_commands_fast(cmds,demultiplex_q,Parallel);
end
    
%-------------------------------------------------------------------------


fprintf(1,'Filter reads... \n') ; tic ;

fparams = {};
cmds = {} ;
dirs ={};
for i=1:length(SampleTable)
    s=SampleTable(i) ;
    cd(s.Sample) ;
    [~,ai] = ismember(s.Alignments,{AlignmentTable.Alignment}) ;
    if any(ai==0)
        error(['No alignment name found in alignment_params.csv: ' s(ai==0).Alignments])        
    end
    [~,fi] = ismember({AlignmentTable(ai).Filter},{FilterTable.Filter}) ;
    if any(fi==0)
        error(['No filter name found in read_filter_params.csv: ' AlignmentTable(ai(fi==0)).Filter])        
    end
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
                else
                    fparams{end+1} =  {[pwd '/' s.Sample '.fastq'], [pwd '/' FilterTable(f).Filter '/filter_reads.fastq'], Phred_offset, FilterTable(f).Method, FilterTable(f).Params};
                end
            end
        end
    end
    cd ..
end

if paired
    run_parallel_matlab_commands('filter_reads_paired',fparams,filter_q,Parallel);
else
    run_parallel_matlab_commands('filter_reads', fparams, filter_q, Parallel);
end

run_parallel_unix_commands_fast(cmds,filter_q,Parallel, dirs);

%-------------------------------------------------------------------------



fprintf(1,'Align... \n') ; tic ;

cmds = {} ;
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
        
        
        %Check to see if reference genome is formatted properly for bowtie2
        if ~exist(['../../Reference_Genomes/' ae.Genome '/genome_bowtie2.1.bt2'],'file')
            fprintf(1,'Format reference genome for bowtie... \n') ; tic ;
            c={}; c{end+1}=['/opt/bowtie2/bowtie2-build genome.fasta genome_bowtie2'];
            d=[]; d{end+1}=['../../Reference_Genomes/' ae.Genome];
            run_parallel_unix_commands_fast(c,'empty',0,d);
        end
        
        
        
        if ~exist([dr '/aligned.sorted.bam'])      
            dirs{end+1} = [pwd '/' dr] ;
            switch ae.Method
                case 'bowtie'
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie/bowtie ' ae.Param1 '--phred64-quals --max multialigned.fastq --un unaligned.fastq ' ...
                            '../../../../Reference_Genomes/' ae.Genome '/genome ../filter_reads.fastq aligned.sam' ] ;
                    else
                        cmds{end+1} = ['/opt/bowtie/bowtie ' ae.Param1 ' --max multialigned.fastq --un unaligned.fastq ' ...
                            '../../../../Reference_Genomes/' ae.Genome '/genome ../filter_reads.fastq aligned.sam' ] ;
                    end

                case 'bowtie2'
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 --phred64 -x ../../../../Reference_Genomes/' ae.Genome '/genome_bowtie2' ...
                            ' -U  ../filter_reads.fastq -S aligned.sam --un unaligned.fastq '];
                    else
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -x ../../../../Reference_Genomes/' ae.Genome '/genome_bowtie2' ...
                            ' -U  ../filter_reads.fastq -S aligned.sam --un unaligned.fastq '];
                    end
                case 'bowtie2paired'
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie2/bowtie2  --phred64  -x ../../../../Reference_Genomes/' ae.Genome '/genome_bowtie2' ...
                            '-1 ../filter_reads_1.fastq -2 ../filter_reads_2.fastq -S aligned.sam'];
                    else
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -x ../../../../Reference_Genomes/' ae.Genome '/genome_bowtie2' ...
                             ' -1 ../filter_reads_1.fastq -2 ../filter_reads_2.fastq -S aligned.sam '];
                    end
                case 'bowtie2pairedfilter' % allows no ambiguous characters
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --un-conc unaligned.fastq --phred64 -x  ../../../../Reference_Genomes/' ae.Genome '/genome_bowtie2' ...
                            '-1 ../filter_reads_1.fastq -2 ../filter_reads_2.fastq -S aligned.sam'];
                    else
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -X 2000 --no-mixed --very-sensitive --n-ceil 0,0.01 --un-conc unaligned.fastq -x  ../../../../Reference_Genomes/' ae.Genome '/genome_bowtie2' ...
                             ' -1 ../filter_reads_1.fastq -2 ../filter_reads_2.fastq -S aligned.sam '];
                    end
                case 'bowtie2pairedxt' % allows no ambiguous characters
                    if Phred_offset == 64
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -X 2000 --no-mixed --dovetail --very-sensitive --n-ceil 0,0.01 --un-conc unaligned.fastq --phred64 -x  ../../../../Reference_Genomes/' ae.Genome '/genome_bowtie2' ...
                            '-1 ../filter_reads_1.fastq -2 ../filter_reads_2.fastq -S aligned.sam'];
                    else
                        cmds{end+1} = ['/opt/bowtie2/bowtie2 -X 2000 --no-mixed --dovetail --very-sensitive --n-ceil 0,0.01 --un-conc unaligned.fastq -x  ../../../../Reference_Genomes/' ae.Genome '/genome_bowtie2' ...
                             ' -1 ../filter_reads_1.fastq -2 ../filter_reads_2.fastq -S aligned.sam '];
                    end
                    
                                        
            end
        end
        all_dirs{end+1} = [pwd '/' dr] ;
        all_genomes{end+1} = ae.Genome ;
    end
    cd ..
end


if ~exist('alignment_stats','dir')
    run_parallel_unix_commands_fast(cmds,alignment_q,Parallel,dirs);
    mkdir('alignment_stats')
    !mv run_parallel_unix_commands_fast_tmp/*  alignment_stats/
end


fprintf(1,'Create .bam and sort... \n') ; 
run_parallel_unix_commands_fast({'/opt/samtools/bin/samtools view -bS -o aligned.bam aligned.sam'},alignment_q,Parallel,dirs);
run_parallel_unix_commands_fast({'/opt/samtools/bin/samtools sort aligned.bam aligned.sorted'},alignment_q,Parallel,dirs);
run_parallel_unix_commands_fast({'rm aligned.sam'},alignment_q,0,dirs);
run_parallel_unix_commands_fast({'rm aligned.bam'},alignment_q,0,dirs);


%-------------------------------------------------------------------------



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
            cmds{end+1} = ['/opt/samtools/bin/samtools mpileup -q30 -6 -O -s -B -d3000 -f  ../../../../Reference_Genomes/' all_genomes{i} '/genome.fasta aligned.sorted.bam > strain.pileup'] ;
        else
            cmds{end+1} = ['/opt/samtools/bin/samtools mpileup -q30 -s -O -B -d3000 -f ../../../../Reference_Genomes/' all_genomes{i} '/genome.fasta aligned.sorted.bam > strain.pileup'] ;
        end
        dirs{end+1} = all_dirs{i} ;
    end
end


run_parallel_unix_commands_fast(cmds,processing_q,Parallel,dirs);





dirs = {} ;
cmds = {} ;
for i=1:length(all_dirs)
   if ~exist([all_dirs{i} '/variant.vcf'],'file') ;
        if Phred_offset==64
            cmds{end+1} = ['/opt/samtools/bin/samtools mpileup -q30 -6 -S -ugf ../../../../Reference_Genomes/' all_genomes{i} '/genome.fasta aligned.sorted.bam > strain'] ;
        else
            cmds{end+1} = ['/opt/samtools/bin/samtools mpileup -q30 -S -ugf ../../../../Reference_Genomes/' all_genomes{i} '/genome.fasta aligned.sorted.bam > strain'] ;
        end
        
        dirs{end+1} = all_dirs{i} ;
    end
end

run_parallel_unix_commands_fast(cmds,processing_q,Parallel,dirs);


%vcf files are called for the purpose of making calls of difference from
%reference
run_parallel_unix_commands_fast({'/opt/samtools/bin/bcftools view -g strain > strain.vcf'},alignment_q,Parallel,dirs);
run_parallel_unix_commands_fast({'/opt/samtools/bin/bcftools view -vS strain.vcf > variant.vcf'},alignment_q,Parallel,dirs);


%-------------------------------------------------------------------------



fprintf(1,'Create diversity.mat... \n') ; tic ;


params = {} ;
for i=1:length(all_dirs)
    if ~exist([all_dirs{i} '/diversity.mat'],'file') ;
        params{end+1} =  {[all_dirs{i} '/strain.pileup'], [all_dirs{i} '/diversity.mat'], all_genomes{i}};
   end
end

run_parallel_matlab_commands('pileup_to_diversity_matrix', params, processing_q, 1);



fprintf(1,'Building .bai files for viewing alignment... \n') ; tic ;


%Rename bam first, then make bai files, then delete temporary files


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

%-------------------------------------------------------------------------



% Untested, might break -- need to adjust integer inputs to coverageMaps.py
% to be appropriate genome size 

% fprintf(1,'Creating coverage maps \n') ; tic ;
% 
% dirs={};
% 
% for i=1:length(all_dirs)
%     if ~exist([all_dirs{i} '/coverage_smoothed_500.png'],'file') ;
%         dirs{end+1}=all_dirs{i};
%     end
% end
% 
%     
% run_parallel_unix_commands_fast({'python2.5 ../../../../scripts/coverageMaps.py 35 22 9'},plotting_q,0,dirs);
% 


