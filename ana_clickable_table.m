function mutgen = ana_clickable_table(StrainNames, Call3, Mut3, MutGenVCF3, ScafNames, RefGenome, ind, z, params, Fig)
% 2012 Feb, Idan Yelin & Roy Kishony

IsGenomeLoaded = false ;

DP4 = reshape([MutGenVCF3.DP4],4,size(MutGenVCF3,1),size(MutGenVCF3,2)) ;

li = length(ind) ;
lz = length(z) ;
lp = length(params) ;

mutgen = {} ;
for i=1:lz
    mutgen{1,i} = StrainNames(z(i)).Sample(end-2:end) ;
end
mutgen(2:li+1,1:lz) = num2cell(Call3(ind,z)) ;

params_s = {
    {'Gene',    {Mut3.gene}     },...
    {'Protein', {Mut3.protein}  },...
    {'Scaf',    {Mut3.scafold}  },...
    {'Pos',     {Mut3.pos}      },...
    {'NonSyn',  {Mut3.NonSyn}   },...
    } ;

params_s = [ params_s, params ] ;

for i=1:length(params_s)
    mutgen{1,lz+i} = params_s{i}{1} ;
    mutgen(2:li+1,lz+i) = params_s{i}{2}(ind) ;
end

% open mutgen
% xlswrite('Mutation_list',mutgen)

%%
figure(Fig);

set(gcf, 'Position',[10         50        1250         550]);
columnname =   mutgen(1,1:end);

%columnformat = {'numeric', 'bank', 'logical', {'Fixed' 'Adjustable'}};
%columneditable =  [false false true true];
t = uitable('Units','normalized','Position',[0 0 1 1], 'Data', mutgen(2:end,:),...
    'ColumnName', columnname,...
    'RowName',[],'ColumnWidth',num2cell([ones(1,lz)*30, 100, 200, 20, 20 ,20, ones(1,lp)*20]), ...
    'CellSelectionCallback',@mut_matix_clicked );
%            'ColumnFormat', columnformat,...
%            'ColumnEditable', columneditable,...


    function mut_matix_clicked(src, event)
        global cMG cM
        % need to install IGV viewer and set it up for communication in Advanced
        % settings. then remover the comments %*
        rc = event.Indices ;
        dt = get(src,'data') ;
        
        t = tcpip('localhost', 60152) ;
        fopen(t) ;
        if rc(2)<=lz
            cMG = MutGenVCF3(ind(rc(1)),z(rc(2))) ;
        end
        cM = Mut3(ind(rc(1))) ;
        open cMG
        open cM
        
        if size(rc,1)>1
            run_cmd('new') ;
        else
            figure(100);clf
            bar(squeeze(DP4(:,ind(rc(1)),z))','stack')
            if rc(2)<=length(z)
                bai_name = '' ;
                bai_name = [StrainNames(z(rc(2))).ExperimentFolder '/' StrainNames(z(rc(2))).Sample '/' StrainNames(z(rc(2))).AlignmentFolder '/aligned.sorted.bam.bai' ]
                    %StrainNames(z(rc(2))).Sample '.bam.bai' ] ;
                
                
                if ~exist(['../' bai_name],'file')
                    error('Create bai files for viewing alignment')
                    % eval(['!/opt/bin/samtools index ' bai_name(1:end-4)])
                end
                
                if ~IsGenomeLoaded
                    run_cmd(['genome  Reference_Genomes/' RefGenome '/' RefGenome '.genome' ])
                    IsGenomeLoaded = true ;
                end
                
                run_cmd(['load ' bai_name(1:end-4)]) ;
                
                run_cmd(['goto ' ScafNames{cM.scafold} ':' num2str(cM.pos)])
                
            end
        end
        
        
        fclose(t);
        delete(t);
        clear t
        
        function run_cmd(c)
            disp(c)
            fprintf(t, c)
            response = fgetl(t)
        end
        
    end


end
