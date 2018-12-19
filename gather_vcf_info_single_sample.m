function gather_vcf_info_single_sample(filename,SampleName)



%Searches for all lines in strain.vcf that match
%that position. Only record info for best position.

%P must be sorted


load('for_gather_vcf_info');

pos_list=chrpos2index(Positions,ChrStarts);

%check to make sure Positions is formated properly
sortedpos_list=sort(pos_list);
if (sum(pos_list==sortedpos_list))~=(numel(pos_list))
    error('list of positions inputed to gather_vcf_info_single_sample needs to be sorted!\n')
elseif (numel(unique(pos_list)))~=(numel(pos_list))
    error('list of positions inputed to gather_vcf_info_single_sample needs to be unique!\n')
end


Tab=9 ;

quals=zeros(numel(pos_list),1);
calls=repmat('N',numel(pos_list),1);



fid = fopen(filename,'r') ;

i=0 ;
ipos=1 ; %position in pos_list we are scanning for
tic
while ~feof(fid) && ~isempty(ipos) %ipos will be empty we are further on genome than last position in search bag
    i=i+1 ;
    if ~mod(i,50000), fprintf(1,'.'); end
    s = fgetl(fid) ;
    f = find(s==Tab,6) ;
    if s(1)~='#' & length(f)>=2 %is a vcf entry?
        
        %get position
        raw_pos = str2double(s(f(1)+1:f(2)-1)) ;
        scf = s(1:f(1)-1) ;
        if isempty(raw_pos)
            raw_pos=0;
        end
        pos=chrpos2index([find(strcmp(scf,ScafNames)), raw_pos] ,ChrStarts);
        %if position is in search bag, store it if the quality is better
        if any(pos==pos_list)
            ipos=find(pos==pos_list,1);
            [c, q] = vcfline2callqual(s);
            if q < quals(ipos)  %more negative is more confidence
                calls(ipos)=c;
                quals(ipos)=q;
            end
            %if position is larger, increment location to check in search bag
        elseif pos>pos_list(ipos)
            ipos = find(pos_list>pos,1) ;
        end
        %otherwise, get the next line
    end
end
fprintf(1,' %6.0f \n',toc) ;


save([TEMPORARYFOLDER '/vcfinfo_' char(SampleName)], 'calls', 'quals');

end

