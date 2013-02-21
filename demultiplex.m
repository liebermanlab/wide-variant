function cmds = demultiplex(IsolateTable)

cmds = {} ;

for i=1:length(IsolateTable)
    s = IsolateTable(i) ;
    if ~exist(s.Sample,'dir')
        mkdir(s.Sample)
    end
    if ~exist([s.Sample '/' s.Sample '.fastq'],'file')        
        bcd_s = s.Barcode ;
        for k=1:length(s.Barcode)
            for n='ACTG'
                if n~=s.Barcode(k)
                    new_barcode = s.Barcode ;
                    new_barcode(k) = n ;
                    bcd_s = [bcd_s, '\|#', new_barcode ] ;
                end
            end
        end        
        cmds{end+1} = [ 'sed -n -e ''/#' bcd_s '/ {N; p}'' ../raw_data/' s.Batch '/s_' num2str(s.Lane) '_sequence.txt > ' s.Sample '/' s.Sample '.fastq' ] ;
    end
end


end