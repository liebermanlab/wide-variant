function chr = scafname2chr(scafname, scafdict, refgenome)
	% HC custom written function for Steno Ab55555
	if nargin < 2
		names = {}; 

		# read genome + extract contig names  
		fastafile = fastaread(['/groups/kishony/Reference_Genomes/' refgenome '/genome.fasta']); 
		headers = {fastafile.Header}; 
		for h = 1:length(headers)
			blanks = strfind(headers{h},' '); % look for blank spaces (for steno Ab55555) 
			if isempty(blanks)
				names{end+1} = headers{h}; 
			else
				names{end+1} = headers{1:blanks(1)}; 
			end 
		end
		numbers = 1:length(headers); 
		
		scafdict = containers.Map(names, numbers); 
	end
	chr = scafdict(scafname); 
end
