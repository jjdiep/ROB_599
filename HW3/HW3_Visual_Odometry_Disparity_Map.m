% load in images
Il = rgb2gray(imread('It1_L.png'))
Ir = rgb2gray(imread('It1_R.png'))
min_d=0;
max_d=64;
w_radius= 3.0; % half of the window size
[D] = genDisparityMap(Il, Ir, min_d, max_d, w_radius)
%

function [D] = genDisparityMap(I1, I2, min_d, max_d, w_radius)
%
% INPUTS:
% I1 = left stereo image
% I2 = right stereo image
% maxdisp = determines the region of interested to search in the right image.
%           It can be thought of as the 'shift' of the right image
% w_radius = half of window size
%
ws = 2*w_radius+1;
% PSEUDO CODE:
% for each pixel (i,j) in the left image:
%     %for each shift value d of the right image:
%     for d = 0:maxdisp
%         % get the template using ws from the left image
%         % compare to the right image, which is shifted by d
%         % compute the Sum of absolute difference, which is a scalar
%         disp_img(i, j, d+1) = SAD
%     end
% end

kernel = ones(ws);
maxdisp = max_d - min_d;

for d = 0:1:maxdisp
    I2_trans = imtranslate(I2,[d 0]);
    I_diff = I1-I2_trans;
    SAD = conv2(I_diff,kernel,'same');
    SAD = abs((SAD));
    disp_img(:,:,d+1) = SAD;
end
% because of our chosen similarity metric, we want to minimize SAD values:
[M,D] = min(disp_img,[],3);%? FILL IN based upon the structure of disp_image

%
% visualize disparity map
figure;
imagesc(D,[0 max_d]);
colormap(gray);
colorbar;
axis image;
%

% EOF
end
