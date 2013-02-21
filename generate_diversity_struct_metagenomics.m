function  [Counts, Windows] = generate_diversity_struct_metagenomics(filename, SampleDirs, p, numfields, starts, totallength)

markersizes=[starts; totallength] - [0; starts]; sizes[1]=[];
maxsize=max(markersizes);

fprintf(1,'Creating counts 3 dimensional matrix \n') ;

Counts=zeros(numfields, length(p), length(SampleDirs));
Windows=nan(maxsize, length(p), length(SampleDirs));

starts(1)=1;

for i=1:length(SampleDirs)
    
    %make 3d matrix counts
    fprintf(1,'Sample: %g  \n',i) ;
    load([SampleDirs{i} '/' filename]);
    Counts(:,:,i)=data(:,p);
    
    %make windows
    for j=1:length(p)
       
        pos=p(j);
        
        %find which marker we're in
        m=find(starts <= pos, 1, 'last');
        mstart=starts(m);
        mend=starts(m+1)-1;
        
        c=data(1:4,mstart:mend)+data(5:8,mstart:mend);
        [sorted, ~] = sort(c,1);
        Windows(mstart:mend,j,i)=(double(sorted(end,:,:))./sum(c,1))';

    end
    

end

