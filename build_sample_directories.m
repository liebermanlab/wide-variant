source_dir = 'G:\PhD\Illumina_Dec2011\Raw_Reads' ;
target_dir = 'G:\PhD\Automated' ;
[~,~,d] = xlsread('isolate_table.xlsx') ;
IsolateTable.Sample = d(:,1) ;
IsolateTable.RefGenome = d(:,2) ;
clear d

for i=1:length(IsolateTable.Sample)
    smp = IsolateTable.Sample{i} ;
    ref = IsolateTable.RefGenome{i} ;
    tdr = [target_dir '\' smp] ;
    sdr = [source_dir '\' smp] ;
    if ~isnan(ref)
        if exist(tdr,'dir')
            fprintf(1,'%g exists. skiping. delete directory to rebuild\n',tdr)
        else
            fprintf(1,'Building directory: %s \n',tdr) ;
            dr=dir([sdr '\*.gz']) ;
            dr = {dr.name} ;
            mkdir(tdr) ;
            fastq_name = {} ;
            bowtie_line = '' ;
            
            fprintf(1,'Unzip, Filter... \n') ;
            for j=1:length(dr)
                cd 'C:\Program Files\7-Zip'
                eval(['!7z e -o' tdr ' ' sdr '\' dr{j}])
                cd(tdr) ;
                fastq_name{j} = filter_reads(dr{j}(1:end-3)) ;
                bowtie_line = [bowtie_line ,',' tdr '\' fastq_name{j}] ;
            end
            
            if ~exist(['Genome_' ref])
                bowtie_line(1) = [] ;
                mkdir(['Genome_' ref])
                fprintf(1,'Bowties... \n') ;
                cd(['Genome_' ref])
                eval([ '!bowtie -S -v 3 -m 1 --max multialigned.fastq --un unaligned.fastq ' target_dir ...
                    '\Reference_Genomes\' ref '\' ref ' ' bowtie_line ' aligned.sam' ]) ;
                
                fprintf(1,'SamTools... \n') ;
                !samtools view -bS -o aligned.bam aligned.sam
                !samtools sort aligned.bam aligned.sorted
                
                fprintf(1,'Pileup... \n') ;
                cyghome = 'c:\cygwin\home\kishonylab\' ;
                copyfile([target_dir '\Reference_Genomes\' ref '\' ref '.txt'], [cyghome 'genome.txt']) ;
                copyfile('aligned.sorted.bam', cyghome) ;
                
                system('C:\cygwin\bin\bash --login -c "samtools mpileup -ugf genome.txt aligned.sorted.bam > strain"') ;
                system('C:\cygwin\bin\bash --login -c "bcftools view -g strain > strain.vcf"') ;
                system('C:\cygwin\bin\bash --login -c "bcftools view -vS strain.vcf > variant.vcf"') ;
                
                movefile([cyghome 'strain.vcf'], '.') ;
                movefile([cyghome 'variant.vcf'], '.') ;
                fprintf(1,'Clean up... \n') ;
                delete([cyghome 'strain'])
                delete([cyghome 'genome.txt'])
                delete([cyghome 'genome.txt.fai'])
                delete([cyghome 'aligned.sorted.bam'])
                cd(tdr) ;
                for j=1:length(dr)
                     delete(dr{j}(1:end-3))
                     delete ([dr{j}(1:end-9) '_filtered.fastq'])
%                      delete(fastq_name{j})
                end
            end
        end
    end
end




