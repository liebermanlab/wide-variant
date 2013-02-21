function [Call,Gen] = build_call(MutGenVCF) 

Call=char(zeros(size(MutGenVCF))+'-') ;
Gen=zeros(size(MutGenVCF)) ;

for i=1:size(MutGenVCF,2)
    for k=1:size(MutGenVCF,1) ;
        alt = MutGenVCF(k,i).alt ;
        ref = MutGenVCF(k,i).ref ;
        if isempty(alt) % No line VCF
            Gen(k,i) = nan ;
            Call(k,i) = 'N' ;
        elseif isnan(alt) % alt='.'
            Gen(k,i) = 0 ;
            Call(k,i) = ref ;
        else
            Gen(k,i) = 1 ;
            if any(alt==',')
                Call(k,i) = 'n' ;
            elseif length(alt)==length(ref)
                Call(k,i) = alt ;
            elseif length(alt)>length(ref)
                Call(k,i) = 'I' ;
            else
                Call(k,i) = 'D' ;
            end
        end
    end
end

return