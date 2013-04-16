function [annotations, allannotations, sorted_tabledata] = div_clickable_table(muts, calls, allp, ancnti, cnts, fwindows, cwindows, mutAF, isdiverse, MutQual, RefGenome, ScafNames, SampleInfo, ChrStarts, promoterdistance, showlegends, QualSort)


%p input is params

fprintf(1,'Generating table...\n')

% 2012 Feb, Tami Lieberman, Idan Yelin & Roy Kishony
scrsz = get(0,'ScreenSize');
strand=['g-','k-'];
window_size=floor(size(fwindows,1)/2);
IsGenomeLoaded = false ;
NTs='ATCG';
show_alignment=1;



%Calculate important things
d=mutAF>0;
[maf, majorNT, minorNT] = div_major_allele_freq(cnts);

goodmaf=zeros(size(d)); goodmaf(d>0)=maf(d>0);



Npositions=numel(allp);
Nsamples=size(maf,2);





%Modify annotations such that it is a stand-alone record
annotations=muts;
for i=1:Npositions
    
    annotations(i).qual=MutQual(i);
    annotations(i).AApos=floor(((double(annotations(i).nt_pos)-1)/3+1)*10)/10;
    
    
    if floor(annotations(i).gene_num)==annotations(i).gene_num
        %intragenic
        protein='';
        for j=1:size(annotations(i).protein,1)
            protein=[protein '' annotations(i).protein(j,:)];
        end
        annotations(i).annotation=protein;
        if numel(annotations(i).AA)~=4
            annotations(i).type='U';
        end
        
    else
        %intergenic
        
        annotations(i).annotation=[num2str(annotations(i).distance1) ...
            '-' annotations(i).locustag1 ' , ' num2str(annotations(i).distance2) ...
            '-' annotations(i).locustag2 ];
        if (~isempty(annotations(i).distance1) && annotations(i).distance1 > -1*promoterdistance && annotations(i).distance1 < 0) || ...
                (~isempty(annotations(i).distance2) && annotations(i).distance2 > -1*promoterdistance && annotations(i).distance2 < 0)
            annotations(i).type='P';
        else
            annotations(i).type='I';
        end
    end
    
end


goodpositions=zeros(Npositions,1);
%Nonsynonymous vs Syonymous, major and minor allele
for i=1:Npositions
    
    annotations(i).maf=nan(Nsamples,1);
    annotations(i).mutAF=nan(Nsamples,1);
    if ancnti(i) > 0
        annotations(i).nts=[NTs(ancnti(i))];
    else
        annotations(i).nts=char([]);
    end
    
    annotations(i).AAs='';
    
    annotations(i).mutAF(1:Nsamples)=mutAF(i,:);
    
    annotations(i).maf(1:Nsamples)=maf(i,:);
    annotations(i).maf(calls(i,:)=='I')=-1;
    annotations(i).maf(calls(i,:)=='D')=-2;
    
    for j=1:Nsamples
        if d(i,j)>0
            goodpositions(i)=1;
            annotations(i).nts(end+1)=NTs(majorNT(i,j));
            if isdiverse(i,j)>0
                annotations(i).nts(end+1)=NTs(minorNT(i,j));
            end
            if numel(annotations(i).AA)==4
                annotations(i).AAs(end+1)=annotations(i).AA(majorNT(i,j));
                if (isdiverse(i,j)>0) && (~strcmp(annotations(i).AA(majorNT(i,j)),annotations(i).AA(minorNT(i,j))))
                    annotations(i).AAs(end+1)=annotations(i).AA(minorNT(i,j));
                    annotations(i).type='N'; %important to do at this step because could be nonsynonymous variation between, but not within, samples
                elseif ancnti(i)> 0 & ~strcmp(annotations(i).AA(majorNT(i,j)),annotations(i).AA(ancnti(i)));
                     annotations(i).AAs(end+1)=annotations(i).AA(majorNT(i,j));
                     annotations(i).type='N'; %important to do at this step because could be nonsynonymous variation between, but not within, samples
               end
            end
        end
    end
    
    
    if numel(annotations(i).AA)==4 && ~strcmp(annotations(i).type,'N')
        annotations(i).type='S';
    end

    %remove duplicates
    annotations(i).nts=unique(annotations(i).nts);
    annotations(i).AAs=unique(annotations(i).AAs);
    
    if numel(annotations(i).AA) > 0
        %make mutations
        annotations(i).muts={};
        rAA=annotations(i).AA(NTs==annotations(i).ref);
        nrAA=annotations(i).AAs(find(annotations(i).AAs~=rAA));
        for j=1:numel(nrAA)
            annotations(i).muts{j}=[rAA num2str(floor(annotations(i).AApos)) nrAA(j)];
        end
    end
    
    if sum(calls(i,:)=='I' & d(i,:),2) > 0
        annotations(i).nts= [annotations(i).nts 'I'];
        annotations(i).type='D';
    elseif sum(calls(i,:)=='D' & d(i,:),2) > 0
        annotations(i).nts= [annotations(i).nts 'D'];
        annotations(i).type='D';
    end
    
end

%for returning
allannotations=annotations;

%actual table
annotations=annotations(goodpositions>0);
allp=allp(goodpositions>0);
goodmaf(~goodpositions,:)=[];
cnts(:,~goodpositions,:)=[];

if numel(fwindows)>1
    fwindows(:,~goodpositions,:)=[];
    cwindows(:,~goodpositions,:)=[];
end

%generate table
if numel(ScafNames)>1
    colnames={'Qual', 'Type','Chr','Pos', 'Locustag', 'Gene','Annotation', 'AApos', 'NTs', 'AAs'};
    widths=num2cell([48, 16, 12, 58, 50, 38, 250, 42, 40, 40, ones(1,numel(SampleInfo))*26]);
else
    colnames={'Qual', 'Type','Pos', 'Locustag','Gene','Annotation', 'AApos', 'NTs', 'AAs'};
    widths=num2cell([48, 16, 58, 50, 38, 250, 45, 40, 40, ones(1,numel(SampleInfo))*26]);
end

nonsamplecols=numel(colnames);
for i=1:Nsamples
    colnames{end+1}=SampleInfo(i).Sample;
end

tabledata=cell(sum(goodpositions),numel(colnames));

for i=1:numel(annotations)
    if numel(annotations(i).locustag)>0
        locustag=annotations(i).locustag(end-4:end);
    else
        locustag=0;
    end
    %generate table
    if numel(ScafNames)>1
        tabledata(i,1:nonsamplecols)=[{[annotations(i).qual]} {annotations(i).type} {annotations(i).scafold} {annotations(i).pos} ...
        {locustag} {annotations(i).gene} {annotations(i).annotation} {annotations(i).AApos} {[annotations(i).nts]} ...
        {[annotations(i).AAs]}];
    else
        tabledata(i,1:nonsamplecols)=[ {[annotations(i).qual]} {annotations(i).type} {annotations(i).pos} ...
        {locustag} {annotations(i).gene} {annotations(i).annotation} {annotations(i).AApos} {[annotations(i).nts]} ...
        {[annotations(i).AAs]} ];
    end
    
    for j=1:Nsamples
        if (annotations(i).mutAF(j) > 0) && (annotations(i).mutAF(j) < 1)
            n=[num2str(annotations(i).mutAF(j)) '0' '0'];
            tabledata(i,nonsamplecols+j)=cellstr(n(2:4));
        elseif (annotations(i).mutAF(j) == 1)
            tabledata(i,nonsamplecols+j)={'1.0'};
        elseif (annotations(i).mutAF(j) == -1)
            tabledata(i,nonsamplecols+j)={'I'};
        elseif (annotations(i).mutAF(j) == -2)
            tabledata(i,nonsamplecols+j)={'D'};
        else
            tabledata(i,nonsamplecols+j)={'0'};
        end
    end
end

% sort table by descending quality val
if QualSort==1
    [sorted_tabledata, sortedpositions] = sortrows(tabledata,1);
else
    sorted_tabledata = tabledata;
    sortedpositions=1:size(tabledata,2);
end
    
%display table
figure();clf;hold on;
set(gcf, 'Position',[10         50        1250         550]);
uicontrol('Style','text','Position',[400 45 120 20],'String','Vertical Exaggeration')
t = uitable('Units','normalized','Position',[0 0 1 .97], 'Data', flipdim(sorted_tabledata,1),...
    'ColumnName', colnames,...
    'RowName',[], ...
    'CellSelectionCallback',@mut_matix_clicked, ...
    'ColumnWidth',widths);
h.checkbox1 = uicontrol('Units','normalized','Style','checkbox','String','Show Alignment in IGV when clicked (must have IGV viewer open already)','Min',0,'Max',1, 'Value',0, 'Position',[0 .97 1 .03]);


    function mut_matix_clicked(src, event)
        
        
        
        
        scrsz = get(0,'ScreenSize');
        strand=['g-','k-'];
        window_size=floor(size(fwindows,1)/2);
        
        rc = event.Indices ;
        dt = get(src,'data') ;
        
        ind=sortedpositions(rc(1));
        
               
        chr=annotations(ind).scafold;
        position= annotations(ind).pos;
                
        disp(annotations(ind))
        
        if rc(2) > nonsamplecols
            sample=rc(2)-nonsamplecols;
            disp(sample)
        else
            sample=1;
        end
        

        
        
        %Bar charts of counts
        div_bar_charts(squeeze(cnts(:,ind,:)), sample, {SampleInfo.Sample})
        
        
        
        %Plot MAF in region neighboring locus
        region=(find(allp>allp(ind)-window_size,1):find(allp<allp(ind)+window_size,1,'last'));

        if numel(fwindows)>1
           div_maf_window(annotations(ind), allp(ind)-ChrStarts(chr), window_size, [], squeeze(fwindows(:,ind,:)), allp(region)-ChrStarts(chr), goodmaf(region,:), {SampleInfo.Sample},sample,showlegends)
           div_cov_window(annotations(ind), allp(ind)-ChrStarts(chr), window_size, squeeze(cwindows(:,ind,:)), {SampleInfo.Sample},sample, showlegends)
        end
        

        
        show_alignment=get(h.checkbox1, 'Value');
        if show_alignment==1
            
            t = tcpip('localhost', 60152) ;
            fopen(t) ;
            if rc(2) > nonsamplecols
%                 bai_name = '' ;
%                 bai_name = ['../../Illumina_pipeline/' StrainNames(z(rc(2))).ExperimentFolder '/' StrainNames(z(rc(2))).Sample '/' StrainNames(z(rc(2))).AlignmentFolder '/aligned.sorted.bam.bai' ]
                bai_name = ['/Volumes/sysbio/KISHONY LAB/illumina_pipeline/' SampleInfo(sample).ExperimentFolder '/' SampleInfo(sample).Sample '/' SampleInfo(sample).AlignmentFolder '/aligned.sorted.bam.bai' ];
                %bai_name = ['../../Tami/Diversity_Illumina_pipeline/' StrainNames(z(rc(2))).ExperimentFolder '/' StrainNames(z(rc(2))).Sample '/' StrainNames(z(rc(2))).AlignmentFolder '/aligned.sorted.bam.bai' ]
                %StrainNames(z(rc(2))).Sample '.bam.bai' ] ;
                
                
                if ~exist(bai_name,'file')
                    error('Create bai files for viewing alignment')
%                    eval(['!/opt/bin/samtools index ' bai_name([1:end-4])])
                end
                
                %exist(['../Reference_Genomes/' RefGenome '/' RefGenome '.genome' ],'file')
                
                if ~IsGenomeLoaded
                    run_cmd(['genome  Reference_Genomes/' RefGenome '/genome.fasta' ])
                    IsGenomeLoaded = true ;
                end
                
                run_cmd(['load "' bai_name(1:end-4) '"']) ;
                
                run_cmd(['goto "' ScafNames{chr} ':' num2str(position) '"'])
                
            end
            
            
            fclose(t);
            delete(t);
            clear t
            
        end
        
        function run_cmd(c)
            disp(c)
            fprintf(t, c)
            response = fgetl(t)
        end
        
        
        
    end




end
