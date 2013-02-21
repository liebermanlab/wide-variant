
%February 2012. 
%
% edited Seungsoo Kim December 2012
% -
%
%Roy Kishony, Tami Lieberman, and Idan Yelin
%Takes unaligned fastq files from Illumina sequencing and produces a
%single table containing all called SNPs 



%Make sure to that dependencies in run_parallel_matlab_commands are
%properly specified


%-----Important parameters that shouldn't be hardcoded in ideal version----

%Note that CoverageMap is hardcoded right now
%In pileup, don't count bases next to an indel

Phred_offset = 33 ; %This should be taken from fastq files, 64 in older versions

% B. dolosa
% ChrStarts=[0,3413073,5589202];
% GenomeLength=6420400;

%S. aureus NCTC 8325
ChrStarts=0;
GenomeLength=2821361;

Directions= ['12'] ; %if not paired, ['1']
Parallel = true ;

fast_q='sysbio_15m';
demultiplex_q='sysbio_12h';
filter_q='sysbio_2h';
alignment_q='sysbio_12h';
processing_q='sysbio_1d'; % processing_q='shared_1d' ;
plotting_q='sysbio_int_2h';

path('../scripts/',path);

%is a Bdolosa genome?
isBdolosa=0;

%overwrite
overwrite=0;

%-------------------------------------------------------------------------


IsolateTable = read_isolate_table ;
FilterTable = read_filter_table ;
AlignmentTable = read_alignment_table ;


%Folders set up manually, no demultiplexing necc

%Old system had multiplexed data
%fprintf(1,'Demultiplex... \n') ; tic ; 
%cmds = demultiplex(IsolateTable);
%run_parallel_unix_commands_fast(cmds,demultiplex_q,Parallel);


fprintf(1,'Filter reads... \n') ; tic ;

fparams = {};
cmds = {} ;
dirs ={};
for i=1:length(IsolateTable)
    s=IsolateTable(i) ;
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
        %SK: if paired
        if strcmp(Directions,'12')
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
            for d=Directions
                fname_in=[pwd '/' s.Sample '_' d '.fastq'];
                fname_out=[pwd '/' FilterTable(f).Filter '/filter_reads_' d '.fastq'];
                if ~exist(fname_out,'file') || overwrite
                    if strfind(FilterTable(f).Method,'nofilter') %if no filter, copy file into subdirectory-- not the best way to do this
                        cmds{end+1}=['cp ../' s.Sample '_' d '.fastq filter_reads_' d '.fastq'];
                        dirs{end+1}=[s.Sample '/nofilter'];
                    else
                        fparams{end+1} =  {[pwd '/' s.Sample '_' d '.fastq'], [pwd '/' FilterTable(f).Filter '/filter_reads_' d '.fastq'], Phred_offset, FilterTable(f).Method, FilterTable(f).Params};
                    end
                end
            end
        end
    end
    cd ..
end

if strcmp(Directions,'12')
    run_parallel_matlab_commands('filter_reads_paired',fparams,filter_q,Parallel);
else
    run_parallel_matlab_commands('filter_reads', fparams, filter_q, Parallel);
end
run_parallel_unix_commands_fast(cmds,filter_q,Parallel, dirs);



fprintf(1,'Align... \n') ; tic ;

cmds = {} ;
dirs = {} ;
all_dirs = {} ;
all_genomes = {} ;


for i=1:length(IsolateTable)
    s=IsolateTable(i) ;
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

        if ~exist([dr '/aligned.sam'])
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
% run_parallel_unix_commands_fast({'/opt/samtools/bin/samtools view -bS -o aligned.bam aligned.sam'},alignment_q,Parallel,dirs);
% run_parallel_unix_commands_fast({'/opt/samtools/bin/samtools sort aligned.bam aligned.sorted'},alignment_q,Parallel,dirs);

run_parallel_unix_commands_fast({'rm aligned.sam'},alignment_q,0,dirs);

%run_parallel_unix_commands_fast({'rm aligned.bam'},alignment_q,0,dirs);





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
           %     cmds{end+1} = ['/opt/samtools/bin/samtools mpileup -q30 -s -O -B -d3000 -f ../../../../Reference_Genomes/Bdolosa/genome.fasta aligned.sorted.bam > strain.pileup'] ;

        end
        dirs{end+1} = all_dirs{i} ;
    end
end


%temp removal -- put back in
% run_parallel_unix_commands_fast(cmds,processing_q,Parallel,dirs);





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

% run_parallel_unix_commands_fast(cmds,processing_q,Parallel,dirs);


%vcf files are called for the purpose of making calls of difference from
%reference


% run_parallel_unix_commands_fast({'/opt/samtools/bin/bcftools view -g strain > strain.vcf'},alignment_q,Parallel,dirs);

% run_parallel_unix_commands_fast({'/opt/samtools/bin/bcftools view -vS strain.vcf > variant.vcf'},alignment_q,Parallel,dirs);





fprintf(1,'Create diversity.mat... \n') ; tic ;


params = {} ;
for i=1:length(all_dirs)
  %  if ~exist([all_dirs{i} '/diversity.mat'],'file') ;
        params{end+1} =  {[all_dirs{i} '/strain.pileup'], [all_dirs{i} '/diversity.mat'], ChrStarts, GenomeLength};
 %  end
end


if isBdolosa
    run_parallel_matlab_commands('read_pileup_Bdolosa', params, processing_q, 1);
else
    run_parallel_matlab_commands('read_pileup_5', params, processing_q, 1);
end

%stop


fprintf(1,'Building .bai files for viewing alignment... \n') ; tic ;


%Rename bam first, then make bai files, then delete temporary files


cmds_bai = {} ;
dirs = {} ;

for i=1:length(IsolateTable)
    s=IsolateTable(i) ;
    for a = s.Alignments'
        ai = find(strcmp({AlignmentTable.Alignment},a)) ;  
        ae = AlignmentTable(ai) ;
        fi = find(strcmp({FilterTable.Filter},ae.Filter)) ;
        dr = [s.Sample '/' FilterTable(fi).Filter '/' ae.Alignment] ;
        if ~exist([dr '/aligned.sorted.bam.bai'])
            cmds_bai{end+1} = ['/opt/samtools/bin/samtools index ' s.Sample '.bam'];
            dirs{end+1} = [pwd '/' dr] ;
        end
    end
end

run_parallel_unix_commands_fast({'/opt/samtools/bin/samtools index aligned.sorted.bam'},alignment_q,Parallel,dirs);

run_parallel_unix_commands_fast(cmds_bai,alignment_q,Parallel,dirs);




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
% run_parallel_unix_commands_fast({'python2.5 ../../../../scripts/coverageMaps.py 35 22 9'},plotting_q,0,dirs);



