function mut_freq = get_fixed_mut_freq(Calls, MutQual, qual_0, ancestor)
    goodpositions=MutQual>qual_0; 
    callsgood = Calls(goodpositions,:);  
    ancestorgood = ancestor(goodpositions); 
    mut_freq = zeros(size(callsgood,1), size(callsgood,2)); 
    NTs = 'ATCGN'; 
    % make it 0 or 1, depending on if it matches "reference" clone 
    for i = 1:length(ancestorgood)
        
        if ancestorgood(i) ~= 0
            anc = NTs(ancestorgood(i)); 
            mut_freq(i, strfind(callsgood(i,:), anc)) = 1; 
        else
            % find majority nucleotide
            maxoccurrence = 0; 
            for n = 1:length(NTs)
                occ = length(strfind(callsgood(i,:),NTs(n))); 
                if occ > maxoccurrence
                    maxoccurrence = occ; 
                    maxNT = NTs(n); 
                end
            end
            mut_freq(i, strfind(callsgood(i,:), maxNT)) = 1; 
            
        end
            
    end
end