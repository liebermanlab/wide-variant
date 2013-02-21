function [lins,indx] = get_specific_line(filename,ScafNames,mut_list,output_file)
Tab=9 ;
y=0 ;

fid = fopen(filename,'r') ;
lins = {} ;
indx = [] ;
pos_list = mut_list(:,1)*1e8 + mut_list(:,2) ;
[pos_list,ind] = sort(pos_list) ;

i=0 ;
ipos=1 ;

tic
while ~feof(fid) && ~isempty(ipos)
    i=i+1 ;
    if ~mod(i,50000), fprintf(1,'.'); end
    s = fgetl(fid) ;
    f = find(s==Tab,2) ;
    if length(f)>=2
        raw_pos = str2num(s(f(1)+1:f(2)-1)) ;
        if isempty(raw_pos), raw_pos=0 ; end
        scf = s(1:f(1)-1) ;
        iscf = find(strcmp(scf,ScafNames)) ;
        pos = iscf*1e8 + raw_pos ;
        if pos>=pos_list(ipos)
            if any(pos==pos_list)
                y = y+1 ;
                lins{y} = s ;
                indx(y) = find(pos_list==pos,1) ;
            end
            ipos = find(pos_list>pos,1) ;
        end
    end
end
fprintf(1,' %6.0f \n',toc) ;
indx=ind(indx) ; 
fclose(fid) ;

if nargin==4
    save(output_file,'lins','indx')
end

end

