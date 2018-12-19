function e = entropym(matrix)

e=matrix.*log2(matrix);

e(isnan(e))=0;

e=sum(e,2);
e=e*-1;
end
