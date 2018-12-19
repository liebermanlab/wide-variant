function scafdict = makescafdict(refgenome, namedelimiter)
	if nargin < 2
		namedelimiter = ' '; 
	end

	names = {}; 
	
	% read genome + extract contig names
	fastafile = fastaread(['/groups/kishony/Reference_Genomes/' refgenome '/genome.fasta']); 
	headers = {fastafile.Header}; 
	for h = 1:length(headers)
		cur_header = headers{h}; 
		blanks = strfind(cur_header, namedelimiter); 
		if isempty(blanks) 
			names{end+1} = cur_header; 
		else 
			names{end+1} = cur_header(1:blanks(1)-1); 
		end
	end
	numbers = 1:length(headers); 
	scafdict = containers.Map(names, numbers); 
		
end
