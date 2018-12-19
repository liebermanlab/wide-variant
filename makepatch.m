function [childy1, childy2] = makepatch(index,parentx,parenty1,parenty2)


global parent;
global genotypes;
global compartmentfreqs;
global compartment_colors_scatter;
global runningploty_byx;
global plotted;
global sizescale;
global de_novo_muts_pt;
global multiplelineagespt;
global isindivergentlineage;
global genomic_positions;

freqs=compartmentfreqs(index,:);
comparts=find(freqs>0);


%which mutations are new
if index==(size(genotypes,1)+1)
    newmuts=[];
elseif parent(index)==0
    newmuts=find(genotypes(index,:));
else
    newmuts=find(genotypes(index,:)-genotypes(parent(index),:));
end


%figure out x axes coordiate placeholders
if index==(size(genotypes,1)+1)
    %ancestor
    x=0;
    x2=sum(genotypes(1,:))+1; %needed an extra placeholder
else
    x=sum(genotypes(index,:));
    x2=x;
end


%figure out right height for plotting
parent_mid_y=parenty2-(parenty2-parenty1)/2;
child_height=sum(freqs(comparts))*(sizescale);
if runningploty_byx(x2) <= parenty1
    runningploty_byx(x2) = parenty1;
    if x~=0 & sum(parent==parent(index))==1 %if only child, center it
        runningploty_byx(x2) = parent_mid_y - child_height/2;
    end
end

%decide on origin of genotype
% parentcomparts=find(compartmentfreqs(parent(index),:);
% comparts=find(freqs>0);
% linew=2;



%plot line
if x~=0
    
    nummuts=numel(newmuts);

    
    if multiplelineagespt & parent(index)==0
        plot([x parentx], [runningploty_byx(x2)+child_height/2 parent_mid_y],'LineStyle', '-', 'LineWidth', 2, 'Color', 'k')
    else
        if numel(newmuts)>2
            [~,b]=ismember([de_novo_muts_pt(newmuts).pos],genomic_positions);
            isindivergentlineage(b)=numel(newmuts);
        end
        
        yincrement=(runningploty_byx(x2)+child_height/2-parent_mid_y)/nummuts;
        xincrement=(x-parentx)/nummuts;
        for i=1:numel(newmuts)
            color=rgb('Black'); linew=2; style='-'; 
%            % disp(newmuts)
%             if de_novo_muts_pt(newmuts(i)).type=='N'
%                 style='--'; 
%             elseif de_novo_muts_pt(newmuts(i)).type=='S'
%                 style=':'; 
%             else
%                 style='-'; 
%             end
%             
%             if de_novo_muts_pt(newmuts(i)).isnewspread==1 & numel(newmuts)<3
%                 color=rgb('Red');
%             elseif de_novo_muts_pt(newmuts(i)).isnewspread==-1 &numel(newmuts)<3
%                 color=rgb('Blue');
%             else
%                 color=rgb('Black');
%             end
            plot([parentx+(i*xincrement) parentx+((i-1)*xincrement)], [parent_mid_y+(i*yincrement) parent_mid_y+((i-1)*yincrement)],'LineStyle', style, 'LineWidth', linew, 'Color', color)
        end

            
    end
end




%make sure future patches account for the fact this this line is there
if x > 1 & (x - parentx) > 1  %increase runningploty_byx for lower values for x so that lines don't cross
    f=polyfit([x parentx],[runningploty_byx(x2)+child_height/2 parent_mid_y],1);
    for i=max(1,(parentx-.5)):(x-1)
        if runningploty_byx(i) < max((i)*f(1)+f(2),(i+.5)*f(1)+f(2))+sizescale*.05;
            runningploty_byx(i) = max((i)*f(1)+f(2),(i+.5)*f(1)+f(2))+sizescale*.05;
            if f(1)>0
                runningploty_byx(i)=runningploty_byx(i)+sizescale*.1;
            end
        end
    end
end

    
    
%plot the actual node
childy1=runningploty_byx(x2);
childy2=runningploty_byx(x2)+child_height;
if numel(comparts)==0 & x~=0
    plot([x x+.5], [runningploty_byx(x2) runningploty_byx(x2)],'k-')
else
    bottom=runningploty_byx(x2);
    %make patches
    for i=1:numel(comparts)
        c=comparts(i);
        height=freqs(c)*sizescale;
        patch([x x x+.5 x+.5],  [runningploty_byx(x2) runningploty_byx(x2)+height runningploty_byx(x2)+height runningploty_byx(x2)], 'k',...
            'FaceColor', rgb(compartment_colors_scatter{c}), 'LineStyle', 'none')
        % text(x,runningploty_byx(x2),num2str(plotted))
        runningploty_byx(x2)=runningploty_byx(x2)+height;
    end
    %make border
    plot([x x],[bottom runningploty_byx(x2)],'LineWidth',2,'Color',rgb('Black'))
    plot([x+.5 x+.5],[bottom runningploty_byx(x2)],'LineWidth',2,'Color',rgb('Black'))
    plot([x x+.5],[bottom bottom],'LineWidth',2,'Color',rgb('Black'))
    plot([x x+.5],[runningploty_byx(x2) runningploty_byx(x2)],'LineWidth',2,'Color',rgb('Black'))
    %right side is thicker to mitigate people thinking that the location of
    %the left side of connecting lines matters
end


%disp(plotted)
plotted=plotted+1;

