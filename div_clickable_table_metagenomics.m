function [annotations, allannotations] = div_clickable_table_metagenomics(muts, cnts, windows, p, SampleNames, ChrStarts)


% 2012 Feb, Tami Lieberman, Idan Yelin & Roy Kishony

scrsz = get(0,'ScreenSize');
strand=['g-','k-'];
window_size=floor(size(dwinds,1)/2);
IsGenomeLoaded = false ;
NTs='ATCG';
show_alignment=1;


%Calculate important things
[maf, majorNT, minorNT] = div_major_allele_freq(cnts);
d=div_test_thresholds(cnts,p);
goodpositions=sum(d,2)>0;


Npositions=numel(p);
Nsamples=size(maf,2);


%generate table
colnames={'Type','Marker#','Pos','MarkerName','NTs', 'AAs',};
nonsamplecols=numel(colnames);
for i=1:Nsamples
    colnames{end+1}=SampleNames(i).Sample;
end


annotations=muts;


%Nonsynonymous vs Syonymous, major and minor allele
for i=1:Npositions
    
    annotations(i).nts='';
    annotations(i).AAs='';    
    annotations(i).maf=dmaf(i,:);
    annotations(i).type='S'; %by default
    
    for j=1:Nsamples
        if d(i,j)>0
            annotations(i).nts(end+1)=NTs(dmajorNT(i,j));
            annotations(i).nts(end+1)=NTs(dminorNT(i,j));
            annotations(i).AAs(end+1)=annotations(i).AA(dmajorNT(i,j));
            if ~strcmp(annotations(i).AA(dmajorNT(i,j)),annotations(i).AA(dminorNT(i,j)))
                annotations(i).AAs(end+1)=annotations(i).AA(dminorNT(i,j));
                annotations(i).type='N'; %important to do at this step because could be nonsynonymous variation between, but not within, samples
            end
        end
    end
    
    %remove duplicates
    annotations(i).nts=unique(annotations(i).nts);
    annotations(i).AAs=unique(annotations(i).AAs);
    
end

allannotations=annotations;

%actual table
annotations=annotations(goodpositions>0);
allp=allp(goodpositions>0);
goodmaf=zeros(size(d)); goodmaf(d>0)=dmaf(d>0); goodmaf(~goodpositions,:)=[];
windows(:,~goodpositions,:)=[];
cnts(:,~goodpositions,:)=[];

tabledata=cell(sum(goodpositions),numel(colnames));
for i=1:numel(annotations)
    %generate table
    tabledata(i,1:nonsamplecols)=[{annotations(i).type} {annotations(i).marker} {annotations(i).pos} ...
        {annotations(i).Header} {[annotations(i).nts]} {[annotations(i).AAs]}];
    for j=1:Nsamples
        tabledata(i,nonsamplecols+j)={annotations(i).maf(j)};
    end
end


%display table
figure();
set(gcf, 'Position',[10         50        1250         550]);
t = uitable('Units','normalized','Position',[0 0 1 1], 'Data', tabledata,...
    'ColumnName', colnames,...
    'RowName',[], ...
    'CellSelectionCallback',@mut_matix_clicked );




    function mut_matix_clicked(src, event)
        
        
         
         scrsz = get(0,'ScreenSize');
         strand=['g-','k-'];
         window_size=floor(size(dwinds,1)/2);
         
         rc = event.Indices ;
         dt = get(src,'data') ;
         
         ind=rc(1);
         chr=annotations(ind).marker;
         position= annotations(ind).pos;
%         
%         
%         disp(ind)
%         disp(annotations(ind))
%         
%         
%         if rc(2) > nonsamplecols
%             sample=rc(2)-nonsamplecols;
%         else
%             sample=6;
%         end
%         
%         
%         %Bar charts of counts
%         div_bar_charts(squeeze(cnts(:,ind,:)), sample, {SampleNames.Sample})
%         
%         
%         
%         %Plot MAF in region neighboring locus
%         region=(find(allp>allp(ind)-window_size,1):find(allp<allp(ind)+window_size,1,'last'));
%         div_maf_window(annotations(ind), window_size, iwinds{ind}, squeeze(dwinds(:,ind,:)), allp(region)-ChrStarts(chr), goodmaf(region,:), {SampleNames.Sample})
%         
%         
%         
%         
%         
%         
%         show_alignment=0;
%         if show_alignment
%             chr=annotations(ind).scafold;
%             position= annotations(ind).pos;
%             
%             
%             t = tcpip('localhost', 60152) ;
%             fopen(t) ;
%             
%             bai_name = [SampleNames(sample).ExperimentFolder '/' SampleNames(sample).Sample '/' SampleNames(sample).AlignmentFolder '/aligned.sorted.bam.bai' ];
%             
%             
%             if ~exist(['../' bai_name ], 'file')
%                 error('Create bai files for viewing alignment')
%                 eval(['!/opt/bin/samtools index ' bai_name(1:end-4)])
%             end
%             
%             if ~IsGenomeLoaded
%                 run_cmd(['genome  Reference_Genomes/' RefGenome '/' RefGenome '.genome' ])
%                 
%                 IsGenomeLoaded = true ;
%             end
%             
%             run_cmd(['load ' bai_name(1:end-4)]) ;
%             run_cmd(['goto ' ScafNames{chr} ':' num2str(position)])
%             
%             
%             fclose(t);
%             delete(t);
%             clear t
%             
%         end
%         
%         function run_cmd(c)
%             disp(c)
%           %  fprintf(t, c)
%            % response = fgetl(t)
%         end
%         
%         
%         
%         
%     end
% 
% 


end
