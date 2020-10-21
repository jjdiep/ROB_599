codeword = get_codeword(randn(1000, 64), 5);

function codeword = get_codeword(feature_all,n_c)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
rng(0);
[idx,codeword] = kmeans(feature_all,n_c);

end