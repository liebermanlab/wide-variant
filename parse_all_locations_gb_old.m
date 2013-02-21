function gb = parse_all_locations_gb(old_gb, Sequence) 

gb = old_gb ;

for i=1:length(gb)
    f = parse_location_gb(gb(i).location) ;
    gb(i).loc1 = f.loc1 ; 
    gb(i).loc2 = f.loc2 ; 
    gb(i).strand = ~f.strand ; 
    if f.strand
        gb(i).Sequence = Sequence(f.loc1:f.loc2) ;
    else
        gb(i).Sequence = seqrcomplement(Sequence(f.loc1:f.loc2)) ;
    end    
end

return 
