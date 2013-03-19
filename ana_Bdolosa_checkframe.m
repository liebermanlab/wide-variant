function f=ana_Bdolosa_checkframe(genename, position, start)

f=0;

if strcmp(genename,'BDAG_00258')
    f=2;  % this adjustment makes the start of the protein out of f, but mutations are later and start of protein is probably wrong
end
if strcmp(genename,'BDAG_00389') 
    f=3; %beginning is out of f, but mutation is later
end
if strcmp(genename,'BDAG_00455') 
    f=1;        %end is out of f, but protein probably starts later (blast), but mutation is earlier
end
if strcmp(genename,'BDAG_00703') 
    f=1; %BLAST shows two reasonable fs, probably a different protein from aa1-300 and one after. mutations found are in later part, so use second f.
end
if strcmp(genename,'BDAG_01028') 
    f=1; %BLAST shows two reasonable domains, with the N terminal domain in the 1st f, and C terminal domain in the 2nd f. Mutations found are in the earlier part of the f
end
if strcmp(genename,'BDAG_01528') 
    f=1; %end is out of f, but mutations are earlier
end
if strcmp(genename,'BDAG_01817') 
    f=2; %beginning is out of f, but mutation is beyond this and second f is supported by blast
end
if strcmp(genename,'BDAG_02722')
    f = 1; %mutation is in beginning, end is out of f
end
if strcmp(genename,'BDAG_02789')
    f = 1; %mutations are in middle (941), edge of protein are out of f
end
if strcmp(genename,'BDAG_02801')
    f = 1; %mutations are  in center, beginning out of f
end
if strcmp(genename,'BDAG_02943')
    %the beginning is in a different f than the end, and we
    %have mutations in both.
    %this is on the forward strand
    %these boundaries are approximate
    if (position - start + 1) < 350
        f = 1; %mutations are  in center, beginning out of f
    else
        f = 3;
    end
end
if strcmp(genename,'BDAG_03149')
    f = 2; %mutations is at end, beginning out of f
end
if strcmp(genename,'BDAG_03302')
    f = 3; %found mutation is at end, most of protein is out of f
end
if strcmp(genename,'BDAG_03456')
    f = 1; %found mutation is at middle (739), end is out of f
end
if strcmp(genename,'BDAG_03460')
    f = 3; %found mutation is at beginning (145), edges are out of f
end
if strcmp(genename,'BDAG_03643')
    f = 3; %found mutation is in the middle (24967!) of a very large protein, most of protein is is f 1
end
if strcmp(genename,'BDAG_03828')
    f = 1; %found mutation is at middle (277), end is out of f 
end
if strcmp(genename,'BDAG_04067')
    f = 3; %found mutation is at end of short peptide (84), front is out of f 
end
if strcmp(genename,'BDAG_04305')
    f = 2; %found mutation is in middle (290), front is out of f
end
if strcmp(genename,'BDAG_04545')
    f = 2; %found mutation is at end (378), front is out of f
end
if strcmp(genename,'BDAG_04603')
    f = 1; %found mutations are at either end (218 & 919), everything is in f
end   
if strcmp(genename,'BDAG_01865')
    f = 2; %found mutations is at start (32) , very start is out of f
end
if strcmp(genename,'BDAG_02645')
    f = 1; %found mutation is in middle (146), end is out of f
end
if strcmp(genename,'BDAG_04309')
    f = 1; %found mutation is in middle (254), everything is in f
end