function dNdS = div_calc_dnds(N,S, expected)


dNdS=log(N./S)./log(expected);


end
