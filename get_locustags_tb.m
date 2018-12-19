function [ltagnumbers, intragenic] =get_locustags_tb(p,codingsequences)



cgenes=genomic_position(codingsequences{1} ,p);

%remove intragenic regions
intragenic=cgenes~=floor(cgenes);
cgenes(intragenic)=[];

ltags={codingsequences{1}(cgenes).locustag}';
ltagnumbers=zeros(size(ltags));

num_starts=cellfun(@locustag2num,ltags);

for i=1:numel(ltags)
    tag=ltags{i};
    if tag(end)>=48 & tag(end)<58 %if the tag ends in a number
        ltagnumbers(i)=str2double(tag(num_starts(i):end)) + 6000*~isempty(strfind(tag,'nc'));
        %If than tag contains 'nc' that this is a noncoding element
        %and added a dumby number to it
    else %the tag could end in a 'c' like it does in TB'
        ltagnumbers(i)=str2double(tag(num_starts(i):end-1))+6000*~isempty(strfind(tag,'nc'));
    end
%     c=find(tag=='c' | tag=='A');
%     
%     
%     if isempty(c)
%         ltagnumbers(i)=str2double(tag(3:end));
%     else
%         ltagnumbers(i)=str2double(tag(3:c-1));
%     end
end

    
end

