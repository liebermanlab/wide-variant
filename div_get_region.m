function region = div_get_region(pos,index, wl)

%April 2012. Returns indices of p that are within wl of index

region=pos(find(pos>pos(index)-wl,1):find(pos<pos(index)+wl,1));

end
