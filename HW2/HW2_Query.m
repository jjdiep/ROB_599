function [d,idx] = query()

kitti = load('kitti_bow.mat');
hist_train = kitti.hist_train;
kdtree = kitti.kdtree;
n_c = kitti.n_c;

img_q = rgb2gray(imread('query.png'));
features = SURF(img_q);
hq = get_hist(kdtree,features,n_c);
d = zeros(length(hist_train),1);
for i = 1:1:length(hist_train)
    d(i) = chi_sq_dist(hq,hist_train(i,:));
end
[dmin,idx] = min(d);
end