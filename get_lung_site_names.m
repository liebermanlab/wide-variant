function unique_sitenames = et_lung_site_names(SampleNames)
    sitenames = cell(length(SampleNames),1); 
    for i = 1:length(SampleNames)
        siten = SampleNames{i}(1 : (strfind(SampleNames{i},'-')-1));
        sitenames{i} = siten;
    end
    
    unique_sitenames = unique(sitenames, 'stable'); 
end