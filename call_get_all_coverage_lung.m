function call_get_all_coverage_lung(samplenames) 
    savedir = '~/Dropbox/projects/kishony-lung/figures/BATCH_1/coverage/';
    groupkey = '-'; 
    gnames = get_unique_sample_groups(samplenames); 
    for i = 1:length(gnames)
        gn = gnames{i}; 
        get_all_coverage([gn '-*']); 
        h = gcf; 
        saveas(h, [savedir '20130929_' gn], 'png'); 
        clf(h); 
    end
    
    function unique_sample_groups = get_unique_sample_groups(samplenames)
        groupnames = {}; 
        for s = 1:length(samplenames)
            sn = samplenames{s}; 
            split_p = strfind(sn, groupkey); 
            groupnames{end+1} = sn(1:split_p-1); 
        end
        unique_sample_groups = unique(groupnames); 
    end
end