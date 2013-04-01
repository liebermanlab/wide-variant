function summarize_alignments(alignment_type)


if strcmp(alignment_type,'bowtie2')
    

    files=dir('alignment_stats/out*.txt');
    files={files.name};

    alignments={};

    for i=1:numel(files)
        file=files{i};
        if ~find(file=='.') | ~strcmp(file(1:3),'out')
            error([file ' is not a valid file name for summarize_alignments. Folder alignment_stats should have only tmp and out files'])
        end
        
        fnumber=file(4:(find(file=='.')-1));
        fn=str2double(fnumber);

        %get alignment type
        fid=fopen(['alignment_stats/tmp' fnumber '.sh']);
        garbage = fgetl(fid); %discard first line anyway
        line = fgetl(fid);
        slashes=find(line=='/');
        alignments(i).alignment=line(slashes(end-1)+1:slashes(end)-1);
        alignments(i).sample=line(slashes(end-3)+1:slashes(end-2)-1);    
        fclose(fid);

        %save results    
        fid=fopen(['alignment_stats/out' fnumber '.txt']);
        line = fgetl(fid);
        alignments(i).totalpairs=str2double(line(1:find(line==' ',1)-1));
        line = fgetl(fid);
        line = fgetl(fid);
        alignments(i).unalignedpairs=str2double(line(5:find(line=='2',1)-2));
        line = fgetl(fid);
        alignments(i).uniquepairs=str2double(line(5:find(line=='(',1)-2));
        line = fgetl(fid);
        alignments(i).nonuniquepairs=str2double(line(5:find(line=='(',1)-2));
    end
    
    
    save('../alignment_stats', 'alignments');
    
else 
    error('Only interprets output from bowtie2 alignments')
end