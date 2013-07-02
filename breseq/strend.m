function output = strend(line,pattern)
%STREND Parses a line of text following a given pattern.
%   O = STREND(L,P) returns the string O within string L that immediately
%   follows pattern P, trimmed to the next tab (ASCII = 9)
%
%   Seungsoo Kim
%   4/5/2013

    temp = strfind(line,pattern);
    
    if ~isempty(temp)
        pos1 = temp(1) + length(pattern);
        pos2 = find(char(line(pos1+1:end))==9);
        if ~isempty(pos2)
            output=strtrim(line(pos1:pos1+pos2(1)));
        else
            output=strtrim(line(pos1:end));
        end
    else
        output='';
    end

end