feature_all = randn(1000, 64); % features from all images in the dataset
n_c = 5;
codeword = get_codeword(feature_all, n_c)

feature = randn(20, 64); % features from one image
kdtree = KDTreeSearcher(codeword);
h = get_hist(kdtree, feature, n_c);
%

function h = get_hist(kdtree,feature,n_c)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
IDX = knnsearch(kdtree,feature);
h = hist(IDX,n_c);
h = h/norm(h,1);
end