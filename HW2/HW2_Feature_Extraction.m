image = rgb2gray(imread('query.png'));
feature = SURF(image);

function feature = SURF(image)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

f = detectSURFFeatures(image);
[feature,valid_points] = extractFeatures(image,f);

end