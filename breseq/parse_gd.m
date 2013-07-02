function [mutations] = parse_gd(file)
%PARSE_GD Parse genome diff (gd) files.
%   M = PARSE_GD(F) outputs a struct M with the mutations in gd file (with
%   annotations) F.
%   
%   MOB (mobile elements) and CON (gene conversion) are not currently
%   processed
%
%   Seungsoo Kim
%   4/5/2013

% mutation classes
mutclasses = {'SNP' 'SUB' 'DEL' 'INS' 'MOB' 'AMP' 'CON' 'INV'};

% open genome diff (gd) file
fid = fopen(file,'r');

% initialize mutations array
mutations = {};

% read first line
l = fgetl(fid);
while ischar(l)
    % skip header
    if strncmp(l,'#',1)
    % mutation
    elseif any(strncmp(l,mutclasses,3))
        % mutation number
        i=length(mutations)+1;
        
        % mutation class
        mutations(i).type = mutclasses{strncmp(l,mutclasses,3)};

        % mutation chromosome
        mutations(i).chr = sscanf(l,'%*s %*d %*s %s');

        % mutation position
        mutations(i).position = sscanf(l,'%*s %*d %*s %*s %d');
        
        % gene information
        mutations(i).gene_name = strend(l,'gene_name=');
        mutations(i).gene_position = strend(l,'gene_position=');
        mutations(i).gene_product = strend(l,'gene_product=');
        mutations(i).locus_tag = strend(l,'locus_tag=');        
        
        % other information
        if strncmp(l,'SNP',3) % single nucleotide polymorphism
            mutations(i).alt = char(sscanf(l,'%*s %*d %*s %*s %*d %s'));
            mutations(i).class = strend(l,'snp_type=');
            mutations(i).size = 0;
        elseif strncmp(1,'SUB',3) % substitution
            mutations(i).alt = sscanf(l,'%*s %*d %*s %*s %*d %s');
            mutations(i).size = sscanf(l,'%*s %*d %*s %*s %d');
        elseif strncmp(l,'DEL',3) % deletion
            mutations(i).size = -sscanf(l,'%*s %*d %*s %*s %*d %d');
        elseif strncmp(l,'INS',3) % insertion
            mutations(i).alt = char(sscanf(l,'%*s %*d %*s %*s %*d %s'));
            mutations(i).size = length(mutations(i).alt);
        elseif strncmp(l,'MOB',3) % mobile element; not currently processed
        elseif strncmp(l,'AMP',3) % amplification
            mutations(i).size = sscanf(l,'%*s %*d %*s %*s %*d %d');
            mutations(i).copy = sscanf(l,'%*s %*d %*s %*s %*d %*d %d');
        elseif strncmp(l,'CON',3) % gene conversion; not currently processed
        elseif strncmp(l,'INV',3) % inversion
            mutations(i).size = sscanf(l,'%*s %*d %*s %*s %*d %d');
        end
    end
    % read next line
    l = fgetl(fid);
end

fclose(fid);
end