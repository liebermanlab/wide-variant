function [annotations, allannotations, diversityannotations, isolateannotations] = div_clickable_table_isolates(muts, dpos, ipos, otherpos, calls, cnts, dwinds, cwinds, iwinds, p, coveragethresholds, RefGenome, ScafNames, SampleNames, ChrStarts, promoterdistance)


% 2012 Feb, Tami Lieberman, Idan Yelin & Roy Kishony
scrsz = get(0,'ScreenSize');
strand=['g-','k-'];
window_size=floor(size(dwinds,1)/2);
IsGenomeLoaded = false ;
NTs='ATCG';
show_alignment=1;



%Calculate important things
[imaf, iMinor imajorNT, iminorNT] = strains_major_allele_freq(calls);
[dmaf, dmajorNT, dminorNT] = div_major_allele_freq(cnts);
d=div_test_thresholds(cnts,p, coveragethresholds);
goodmaf=zeros(size(d)); goodmaf(d>0)=dmaf(d>0);
allp=sort(unique([dpos; ipos; otherpos;]));

Npositions=numel(allp);
Nisolates=size(calls,2);
Nsamples=size(dmaf,2);







%Modify annotations in a straightforwards way that doesn't depend on
%calls/counts

% .AA holds the amino acids coresponding to ATCG
% .AAs holds the amino acids actually found

annotations=muts;
for i=1:Npositions
    
    
    annotations(i).AApos=(double(annotations(i).nt_pos)-1)/3+1;
    
    
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
            '-' annotations(i).gene1 ' , ' num2str(annotations(i).distance2) ...
            '-' annotations(i).gene2 ];
        if (~isempty(annotations(i).distance1) && annotations(i).distance1 > -1*promoterdistance && annotations(i).distance1 < 0) || ...
                (~isempty(annotations(i).distance2) && annotations(i).distance2 > -1*promoterdistance && annotations(i).distance2 < 0)
            annotations(i).type='P';
        else
            annotations(i).type='I';
        end
    end
    
end


%Nonsynonymous vs Syonymous, major and minor allele


goodpositions=zeros(Npositions,1);
for i=1:Npositions
    
    annotations(i).maf=nan(Nsamples+1,1);
    annotations(i).nts='';
    annotations(i).AAs='';
    annotations(i).calls=char(zeros(Nisolates,1)) ;
    
    %diverse positions in deep
    annotations(i).maf(2:Nsamples+1)=dmaf(i,:);
    for j=1:Nsamples
        if d(i,j)>0
            goodpositions(i)=1;
            annotations(i).nts(end+1)=NTs(dmajorNT(i,j));
            annotations(i).nts(end+1)=NTs(dminorNT(i,j));
            if numel(annotations(i).AA)==4
                annotations(i).AAs(end+1)=annotations(i).AA(dmajorNT(i,j));
                if ~strcmp(annotations(i).AA(dmajorNT(i,j)),annotations(i).AA(dminorNT(i,j)))
                    annotations(i).AAs(end+1)=annotations(i).AA(dminorNT(i,j));
                    annotations(i).type='N'; %important to do at this step because could be nonsynonymous variation between, but not within, samples
                end
            end
        end
    end
    
    
    %add diverse in isolates
    [fi,iindex]=ismember(allp(i),ipos);
    if fi
        goodpositions(i)=1;

        annotations(i).nts(end+1)=NTs(imajorNT(iindex));
        annotations(i).nts(end+1)=NTs(iminorNT(iindex));
        annotations(i).maf(1)=imaf(iindex);
        annotations(i).calls=calls(iindex,:);
        if numel(annotations(i).AA)==4
            annotations(i).AAs(end+1)=annotations(i).AA(imajorNT(iindex));
            if ~strcmp(annotations(i).AA(imajorNT(iindex)),annotations(i).AA(iminorNT(iindex)))
                annotations(i).AAs(end+1)=annotations(i).AA(iminorNT(iindex));
                annotations(i).type='N';
            end
        end
    end
    
    
    % other positions
    if find(otherpos==allp(i),1) & numel(unique(dmajorNT(i,2:end)))>1 %added this last part in because there are some things that fell through in here
        goodpositions(i)=1;
        annotations(i).nts(end+1:end+6)=NTs(dmajorNT(i,2:end));
        major=mode(dmajorNT(i,2:end)); 
        minor=unique(dmajorNT(i,2:end));
        minor(minor==major)=[];
        if numel(annotations(i).AA)==4
            
            annotations(i).AAs(end+1)=annotations(i).AA(major);

            if numel(minor)>1
                minor=minor(1);
            end
            
             if ~strcmp(annotations(i).AA(major),annotations(i).AA(minor))
                annotations(i).AAs(end+1)=annotations(i).AA(minor); %must specify 1 here because sometimes there are two minor variants
                annotations(i).type='N'; %important to do at this step because could be nonsynonymous variation between, but not within, samples
            end
        end
        
    end
    
    if numel(annotations(i).AA)==4 && ~strcmp(annotations(i).type,'N')
        annotations(i).type='S';
    end
    
    %remove duplicates
    annotations(i).nts=unique(annotations(i).nts);
    annotations(i).AAs=unique(annotations(i).AAs);
       
    %Create 'muts', include synonymous
    if numel(annotations(i).AAs) > 0 %is it a coding position?
        %make mutations
        annotations(i).muts={};
        rAA=annotations(i).AA(NTs==annotations(i).ref);
        othernts=annotations(i).nts; 
        for j=1:numel(othernts)
            if othernts(j)~=annotations(i).ref
                annotations(i).muts{end+1}=[rAA num2str(floor(annotations(i).AApos)) annotations(i).AA(NTs==othernts(j))];
            end
        end
    end
    
end



%for returning
[~,z1]=ismember(dpos,allp);
[~,z2]=ismember(ipos,allp);


diversityannotations=annotations(z1);
isolateannotations=annotations(z2);

allannotations=annotations;




%actual table
annotations=annotations(goodpositions>0);
allp=allp(goodpositions>0);
goodmaf=zeros(size(d)); goodmaf(d>0)=dmaf(d>0); goodmaf(~goodpositions,:)=[];
dwinds(:,~goodpositions,:)=[];
cwinds(:,~goodpositions,:)=[];
iwinds(~goodpositions,:)=[];
cnts(:,~goodpositions,:)=[];




%generate table
colnames={'Type','Chr','Pos','Gene','Annotation', 'AApos', 'NTs', 'AAs', 'IsolatesMAF'};
nonsamplecols=numel(colnames);
for i=1:Nsamples
    colnames{end+1}=SampleNames(i).Sample;
end

tabledata=cell(sum(goodpositions),numel(colnames));

for i=1:numel(annotations)
    %generate table
    tabledata(i,1:nonsamplecols)=[{annotations(i).type} {annotations(i).scafold} {annotations(i).pos} ...
        {annotations(i).gene} {annotations(i).annotation} {annotations(i).AApos} {[annotations(i).nts]} ...
        {[annotations(i).AAs]} {annotations(i).maf(1)}];
    for j=1:Nsamples
        tabledata(i,nonsamplecols+j)={annotations(i).maf(j+1)};
    end
end




%display table
figure();
set(gcf, 'Position',[10         50        1250         550]);
t = uitable('Units','normalized','Position',[0 0 1 1], 'Data', tabledata,...
    'ColumnName', colnames,...
    'RowName',[], ...
    'CellSelectionCallback',@mut_matix_clicked );


% 
% %display table for isolates
% 
% itabledata=cell(sum(goodpositions),numel(colnames));
% for i=1:numel(isolateannotations)
%     %generate table
%     itabledata(i,1:nonsamplecols)=[{isolateannotations(i).type} {isolateannotations(i).scafold} {isolateannotations(i).pos} ...
%         {isolateannotations(i).gene} {isolateannotations(i).annotation} {[isolateannotations(i).nts]} ...
%         {[isolateannotations(i).AAs]} {isolateannotations(i).maf(1)}];
%     for j=1:Nsamples
%         itabledata(i,nonsamplecols+j)={isolateannotations(i).maf(j+1)};
%     end
% end
% figure();
% set(gcf, 'Position',[10         50        1250         550]);
% t = uitable('Units','normalized','Position',[0 0 1 1], 'Data', itabledata,...
%     'ColumnName', colnames,...
%     'RowName',[], ...
%     'CellSelectionCallback',@mut_matix_clicked );




    function mut_matix_clicked(src, event)
        
        
        
        scrsz = get(0,'ScreenSize');
        strand=['g-','k-'];
        window_size=floor(size(dwinds,1)/2);
        
        rc = event.Indices ;
        dt = get(src,'data') ;
        
        ind=rc(1);
        chr=annotations(ind).scafold;
        position= annotations(ind).pos;
        
      %  disp(ind)
        disp(annotations(ind))
        
        
        if rc(2) > nonsamplecols
            sample=rc(2)-nonsamplecols;
        else
            sample=6;
        end
        
        
       % disp(squeeze(cnts(:,ind,:)))
        %Bar charts of counts
        div_bar_charts(squeeze(cnts(:,ind,:)), sample, {SampleNames.Sample})
        
 
        
        %Plot MAF in region neighboring locus
        region=(find(allp>allp(ind)-window_size,1):find(allp<allp(ind)+window_size,1,'last'));
      %  div_maf_window(annotations(ind), allp(ind)-ChrStarts(chr), window_size, iwinds{ind}, squeeze(dwinds(:,ind,:)),  allp(region)-ChrStarts(chr), goodmaf(region,:), {SampleNames.Sample},sample ,0)
       % div_cov_window(annotations(ind), allp(ind)-ChrStarts(chr), window_size, squeeze(cwinds(:,ind,:)), {SampleNames.Sample},sample ,  0)

       samples=[     6     2     5     7     4];
       disp(dmajorNT(ind,samples));
       disp(goodmaf(ind,samples));
        
        
        
        
        show_alignment=0;
        if show_alignment
            chr=annotations(ind).scafold;
            position= annotations(ind).pos;
            
            
            t = tcpip('localhost', 60152) ;
            fopen(t) ;
            
            bai_name = [SampleNames(sample).ExperimentFolder '/' SampleNames(sample).Sample '/' SampleNames(sample).AlignmentFolder '/aligned.sorted.bam.bai' ];
            
            
            if ~exist(['../' bai_name ], 'file')
                error('Create bai files for viewing alignment')
                eval(['!/opt/bin/samtools index ' bai_name(1:end-4)])
            end
            
            if ~IsGenomeLoaded
                run_cmd(['genome  Reference_Genomes/' RefGenome '/' RefGenome '.genome' ])
                
                IsGenomeLoaded = true ;
            end
            
            run_cmd(['load ' bai_name(1:end-4)]) ;
            run_cmd(['goto ' ScafNames{chr} ':' num2str(position)])
            
            
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
