function igv_caller(src, event)
global MutGenVCFind 

target_dir = 'G:\\PhD\\Automated\\' ;


ind = event.Indices ;
dt = get(src,'data') ;
ud = get(src,'userdata') ;
        % ud{1} SNPFileNames ;
        % ud{2} RefGenome
        % ud{3} ScafNames
        % ud{4} Is Genome Loaded?
        % ud{5} DP4

t = tcpip('127.0.0.1', 60151) ;
fopen(t) ;


if size(ind,1)>1
    fprintf(t, 'new') ;
    response = fgetl(t)
else
    if ind(2)==size(dt,2)-6 
        figure(100);clf
        bar(squeeze(ud{5}(:,ind(1),:))','stack')
        
    else
        
        bai_name = [target_dir ud{1}{ind(2)} '\\Genome_' ud{2} '\\aligned.sorted.bam.bai' ] ;
        
        if ~exist(bai_name,'file')
            eval(['!samtools index ' bai_name(1:end-4)])
        end
        
        if ~ud{4}
            disp(['genome ' ud{2} ])
            fprintf(t, ['genome ' ud{2} ]) ;
            response = fgetl(t)
            ud{4} = true ;
        end
        
        disp(['load ' bai_name(1:end-4)])
        fprintf(t, ['load ' bai_name(1:end-4)])
        response = fgetl(t)
        
        
        disp(['goto ' ud{3}{dt{ind(1),end-1}} ':' num2str(dt{ind(1),end})])
        fprintf(t, ['goto ' ud{3}{dt{ind(1),end-1}} ':' num2str(dt{ind(1),end})]  ) ;
        response = fgetl(t)
        
    end
end

fclose(t);
delete(t);
clear t

end

