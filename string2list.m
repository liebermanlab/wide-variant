function values = string2list(s, c)

%s is a string of #s seperated by commas
%c is a one character delimiter


values= sscanf(strrep(s,c,' '),'%f');