
function piehandles = movepieto( piehandles, newx, newy, scalingx, scalingy)
%assume pairs, patch first then text

for K = 1:2:length(piehandles)
    set(piehandles(K), 'Vertices', ...
        bsxfun(@plus, get(piehandles(K), 'Vertices')*[scalingx 0; 0 scalingy], [newx newy]) );
    set(piehandles(K+1), 'Position', get(piehandles(K+1),'Position') + [newx, newy, 0]);
end

end