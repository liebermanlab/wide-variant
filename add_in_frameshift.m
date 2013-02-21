function genes=add_in_frameshift(cds, features)

%Tami Lieberman
%May 2012




%get all genes
allgenes=extract_feature(features,'gene'); 
allgenes=gene_name_from_text(allgenes); %neccessary fix to add gene number to gene field

%convert to gene numbers
allgenesN=div_get_gene_numbers(allgenes);
cdsN=div_get_gene_numbers(cds);

%find overlap
[found,indices]=ismember(allgenesN, cdsN);

genes=allgenes;



for i=1:numel(found);
    if found(i)>0
        genes(i).location=cds(indices(i)).location;
        genes(i).product=cds(indices(i)).product;
        genes(i).codon_start=cds(indices(i)).codon_start;
        genes(i).indices=cds(indices(i)).indices;
        genes(i).protein_id=cds(indices(i)).protein_id;
        genes(i).db_xref=cds(indices(i)).db_xref;
        genes(i).note=cds(indices(i)).note;
        genes(i).translation=cds(indices(i)).translation;
        genes(i).text=cds(indices(i)).text;
        genes(i).Sequence=cds(indices(i)).Sequence;
    end
end





%function below are taken directly from GENBANKREAD

%-------------------------------------------------------------------------%
function theStruct = extract_feature(theText,theFeature)
% Extract the feature information
% typically CDS, gene, mRNA
theCellstr = cellstr(theText);

% look for the tag at the start of a line
feat = strmatch(theFeature,strtrim(theCellstr)); 

featCount = numel(feat);

if featCount == 0
    theStruct = [];
    return;
end

% As the items in the CDS fields is unknown, we rely on indentation to tell
% when the CDS field ends. This is not particularly robust.
indent = find(theCellstr{feat(1)} == theFeature(1),1)-1;
if indent>0
    theCellstr = cellstr(theText(:,indent+1:end));
end
%look for lines with first level features.
featureLines = find(cellfun('isempty',regexp(theCellstr,'^\s')));
featureLines(end+1) = numel(theCellstr)+1;

% Sometimes 'CDS' show up in the /note as a single line
theFeat = intersect(featureLines,feat);
featCount= numel(theFeat);
% create empty struct
theStruct(featCount).location = '';
theStruct(featCount).gene = '';
theStruct(featCount).product = '';
theStruct(featCount).codon_start = [];
theStruct(featCount).indices = [];
theStruct(featCount).protein_id = '';
theStruct(featCount).db_xref = '';
theStruct(featCount).note = '';
theStruct(featCount).translation = '';
theStruct(featCount).text = '';

% loop over all of the CDS
for featloop = 1:featCount

    featurePos = find(featureLines == theFeat(featloop));

    endLine = featureLines(featurePos+1)-1;
    textChunk = strtrim(theCellstr(theFeat(featloop):endLine));

    numLines = numel(textChunk);
    theStruct(featloop).text = char(textChunk);

    textChunk{1} = strtrim(strrep(textChunk{1},theFeature,''));
    [featLocation, startLoop] = getFullText(textChunk, 1);
    theStruct(featloop).location = char(strread([featLocation{:}], '%s'));
    theStruct(featloop).indices = featurelocation(theStruct(featloop).location);

    strLoop = startLoop;
    while(strLoop < numLines)
        [fullstr, endpos] = getFullText(textChunk, strLoop);
        lines = size(fullstr, 1);
        [token,rest] = strtok(fullstr{1},'='); 
        rest= strrep(rest,'"','');
        fullstr{1}= rest(2:end);
        if(lines > 1)
            fullstr{lines} = strrep(fullstr{lines}, '"','');
        end
        switch token(2:end)
            case 'gene'
                theStruct(featloop).gene = char(fullstr{:});
            case 'product'
                theStruct(featloop).product =  char(fullstr{:});
            case 'codon_start'
                theStruct(featloop).codon_start =  char(fullstr{:});
            case 'protein_id'
                theStruct(featloop).protein_id = char(fullstr{:});
            case 'db_xref'
                theStruct(featloop).db_xref = char(fullstr{:});
            case 'note'
                theStruct(featloop).note = char(fullstr{:});
            case 'translation'
                theStruct(featloop).translation = char(strread([fullstr{:}],'%s'));
            otherwise % There may be other fields...
                % disp(sprintf('Unknown field: %s',token));
        end
        strLoop = endpos;
    end
end
%-------------------------------------------------------------------------%
function [fullText, endPos] = getFullText(textChunk, i)
% next qualifier (if exists) starts with '/'
nextKey = find(strncmp('/',textChunk(i+1:end),1));
if isempty(nextKey)
    endPos = numel(textChunk)+1;
else
    endPos = i + nextKey(1);
end
fullText = textChunk(i:endPos-1);



