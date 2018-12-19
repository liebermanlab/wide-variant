x=pwd;
path('~/Dropbox/illumina_pipeline',path);
cd ~/Dropbox/TB' diversity'/Lab' notebook'/_Log' spreadsheets'/;

pts=[730,731,736,737,738,739,742,746,747,755];

donotinclude={'P4-N19-S15','P5-N10-S22', 'P5-N05-S17','P4-N14-S11','P3-N15-S14', 'P5-N22-S03','P5-N07-S20','P5-N03-S24','P5-N02-S24','P5-N02-S21','P4-N20-S10','P5-N11-S23'};

oldnames={};
newnames={};

%% Make plate structure
d=tdfreadunix('Plate_log_simple.csv',',');

k = {} ;
for f=fieldnames(d)'
    if ischar(d.(f{1}))
        d.(f{1}) = cellstr(d.(f{1})) ;
    else
        d.(f{1}) = num2cell(d.(f{1})) ;
    end
    k(end+1,:) = {d.(f{1}){:}} ;
end
PlateTableTemp = cell2struct(k,fieldnames(d),1) ;

for i=1:numel(PlateTableTemp)
    if PlateTableTemp(i).Plate>0
        PlateTable(PlateTableTemp(i).Plate)=PlateTableTemp(i);
    end
end

    
    
%% Make sample structure
d=tdfreadunix('Sample_log_simple.csv',',');

k = {} ;
for f=fieldnames(d)'
    if ischar(d.(f{1}))
        d.(f{1}) = cellstr(d.(f{1})) ;
    else
        d.(f{1}) = num2cell(d.(f{1})) ;
    end
    k(end+1,:) = {d.(f{1}){:}} ;
end
SampleTableTemp = cell2struct(k,fieldnames(d),1) ;



for i=1:numel(SampleTableTemp)
    if SampleTableTemp(i).Pool>0 
        n=SampleTableTemp(i).Column-1+PlateTable(SampleTableTemp(i).Plate).N;
        s=SampleTableTemp(i).Row-65+PlateTable(SampleTableTemp(i).Plate).S;

        %custom changes for plate 6
        if SampleTableTemp(i).Plate==6
            n=24-SampleTableTemp(i).Column+1;
            if SampleTableTemp(i).Column>2
                s=s-1;
            end
        end
        
        if n < 10
            N= ['0' num2str(n)];
        else
            N=num2str(n);
        end
        if s < 10
            S= ['0' num2str(s)];
        else
            S=num2str(s);
        end
        pool=num2str(SampleTableTemp(i).Pool);
        oldnames{i}=['P' pool '-N' N '-S' S];
        newnames{i}=[num2str(SampleTableTemp(i).Patient) '-' SampleTableTemp(i).Sample];
    end
end

nameMap = containers.Map(oldnames, newnames);

save('namemap', 'nameMap', 'oldnames')

c=0;d=0;
 for i=1:length(SampleNames)
     if ismember(SampleInfo(i).Sample, oldnames) & ~ismember(SampleInfo(i).Sample, donotinclude)
        [a,b]=ismember(SampleInfo(i).Sample, oldnames);
        SampleInfo(i).NewName=nameMap(SampleInfo(i).Sample);
        SampleInfo(i).Plate=SampleTableTemp(b).Plate;
        d=d+1;
     else
        disp(SampleInfo(i).Sample)
        c=c+1;
     end
 end
 

sample_ids={};
for i=1:numel(SampleNames)
    if ismember(SampleNames{i},oldnames)  & ~ismember(SampleInfo(i).Sample, donotinclude) & SampleInfo(i).Plate~=8 & SampleInfo(i).Plate~=9;
        sample_ids{i}=nameMap(SampleNames{i});
    else
        sample_ids{i}=SampleNames{i};
    end
end





isolatelist={};
toinclude=zeros(numel(SampleInfo),1);
istrainstudy=zeros(numel(SampleInfo),1);

isolatelist{numel(pts)}=[];


for i=1:numel(pts)
    isolatelist{i}=[];
end

for i=1:numel(SampleInfo)
    s=sample_ids{i};
    if strcmp(s(end-4:end),'study')
        istrainstudy(i)=1;
    elseif ~isempty(str2num(s(1:3))) & ismember(str2num(s(1:3)),pts)
        isolatelist{find(pts==str2num(s(1:3)))}(end+1)=i;
        toinclude(i)=1;
    end
end

save('sample_identity','istrainstudy','isolatelist','toinclude')






 
cd(x);


