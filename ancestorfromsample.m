function [ancnt, ancnti] = ancestorfromsample(calls,sample)


[~,callsi]=ismember(calls,'ATCG');
callsi_mode=mode(callsi,2);

NTs='ATCG';
ancnt=calls(:,sample);
[ancATCG,ancnti]=ismember(ancnt,NTs);
ancnti(~ancATCG)=callsi_mode(~ancATCG);
ancnt(ancnti~=0)=NTs(ancnti(ancnti~=0));
ancnt(ancnti==0)='N';

end