function type=isTrTv(a, b)

%Takes a set of nucleotides (A=1 T=2 C=3 G=4) and returns 1 if the
%mutations is a transversion, 2 if it is a transtion
%If inputs are vectors, each set is made up of element a(i) and b(i) 
%Order of ancestral/derived does not matter


%Method-- hash by multiplying (order does not matter and same is not a possibility!)

   
%transitions
%AG/GA = 4
%CT/TC = 6
    
%transversions
%All others
%AT/TA = 2
%GC/CG = 12
%AC/CA = 3
%GT/TG = 8



type=zeros(size(a));
pair=a.*b;

type(pair>1)=1;
type(pair==4 | pair==6)=2;


%protect against pair being the same

type(a==b)=0;

end
