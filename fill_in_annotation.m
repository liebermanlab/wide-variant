function complete_annotations = fill_in_annotation(annotations)
    python_script_Entrez_fetch = 'python /Volumes/sysbio/kishonylab/illumina_pipeline/scripts-hattie/fetch_gene_annotation.py'; 
    strbit = 'GeneID:'; 
    fieldstosearch = {'annotation', 'protein', 'protein1', 'protein2'}; 
    
    complete_annotations = cell(numel(annotations),1); 
    
    for i = 1:numel(annotations)
        gene = annotations(i).gene;
        if isempty(gene)
            % get Gene ID
            for f = 1:length(fieldstosearch)
                if strfind(annotations(i).(fieldstosearch{f}), strbit)
                    a = annotations(i).(fieldstosearch{f});
                end
            end
            
            if exist('a')
                loc = strfind(a, strbit);
                geneid = a(loc+length(strbit):end); 
                fprintf('\nFetching Gene ID %s\n', geneid); 
                % Fetch Gene ID record from NCBI 
                python_cmd = [python_script_Entrez_fetch ' ' geneid]; 
                [status, gene] = system(python_cmd);
            else
                gene = 'Unknown'; 
            end
            
            clear a
            
        end
        
        % add to complete annotations
        if isempty(annotations(i).muts)
            complete_annotations{i} = [annotations(i).type ' | ' annotations(i).locustag ' | ' gene]; 
        else
            complete_annotations{i} = [annotations(i).type ' | ' annotations(i).muts{1} ' | ' ...
                                        annotations(i).locustag ' | ' gene]; 
        end
    end
end