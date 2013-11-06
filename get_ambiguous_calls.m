function isN = get_ambiguous_calls(Callsgood)
    
    % 1. Compute which positions have too many "N"s
    % 2. Compute which isolates have too many "N"s
    % 3. For each N, get the ratio between the nucleotides called 
    % ___ CASE 1: not enough coverage
    % ___ CASE 2: split evently 
    % ___ CASE 3: should force a call
    
    [num_pos, num_isolates] = size(Callsgood); 
    isN = zeros(num_pos, num_isolates); 
    
    for p = 1:num_pos
        calls_p = Callsgood(p,:);
        Ns = strfind(calls_p,'N'); 
        isN(p,Ns) = 1;
%         num_Ns = length(Ns);
    end
end