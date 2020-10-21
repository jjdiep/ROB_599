% load in the images
It0_left = rgb2gray(imread('It1_L.png'));
It0_right = rgb2gray(imread('It1_R.png'));
It1_left = rgb2gray(imread('It2_L.png'));
It1_right = rgb2gray(imread('It2_R.png'));
%
[R_total,t_total] = motionEstimation3Dto3D(It0_left, It0_right, It1_left, It1_right)
%

function [R,t] = motionEstimation3Dto3D(It0_left, It0_right, It1_left, It1_right)
    %
    % 3D to 3D motion estimation. t0 and t1 refer to the time frames t and
    % t+1, respectively.
    %
    % read calibration parameters (intrinsic and extrinsic) from the
    % calibration file
    calibname = 'calib.txt';
    T = readtable(calibname, 'Delimiter', 'space', 'ReadRowNames', true, 'ReadVariableNames', false);
    B = table2array(T);
    Pleft = vertcat(B(1,1:4), B(1,5:8), B(1,9:12)); % left camera [3x4] projection matrix
    Pright = vertcat(B(2,1:4), B(2,5:8), B(2,9:12)); % right camera [3x4] projection matrix
    %
    % calculate correspondences between It0_left, It0_right stereo pair
%     [MP_L0, MP_R0, MPf_L0, MPf_R0] = detectFeatureMatches(It0_left,It0_right);
%     % calculate correspondences between It1_left, It1_right stereo pair
%     [MP_L1, MP_R1, MPf_L1, MPf_R1] = detectFeatureMatches(It1_left,It1_right);
%     % calculate correspondences between It0_left and It1_left (loc_ft0,loc_ft1)
%     [loc_ft0, loc_ft1, loc_ft0_f, loc_ft1_f] = detectFeatureMatches(It0_left,It1_left);
%     % calculate 3D point cloud Wt0 from It0_left, It0_right stereo pair
%     Wt0 = gen_pointcloud(It0_left, It0_right, Pleft, Pright);
%     % calculate 3D point cloud Wt1 from It1_left, It1_right stereo pair
%     Wt1 = gen_pointcloud(It1_left, It1_right, Pleft, Pright);
    [Wt0, loc_ft0, ft0] = gen_pointcloudandfeatures(It0_left, It0_right, Pleft, Pright);
    [Wt1, loc_ft1, ft1] = gen_pointcloudandfeatures(It1_left, It1_right, Pleft, Pright);

    % Matching by features for Wt0 and Wt1
    indexPairs = matchFeatures(ft0,ft1);
    Wt0 = Wt0(indexPairs(:,1),:);
    loc_ft0 = loc_ft0(indexPairs(:,1),:);
    Wt1 = Wt1(indexPairs(:,2),:);
    loc_ft1 = loc_ft1(indexPairs(:,2),:);
    
    % compute R_total, t_total using ransac
    [R,t] = estimateTransform(Wt0,Wt1,loc_ft0,loc_ft1, Pleft);
%EOF
end 
function [R,t] = estimateTransform(Wt0,Wt1,loc_ft0,loc_ft1, Pleft)
%
%   Find the R,t between source and target point clouds (Wt0 and wt1, respectively)
%   using a set of corresponding points whose locations are given by loc_ft0,loc_ft1. 
%   RANSAC method is used.
%
% OUTPUTS:
%   The estimated R and t
%
% Extract the RANSAC coefficients
coeff.numPtsToSample = 3;     % the minimum number of correspondences needed to fit a model
coeff.iterNum = 2000;         % number of iterations to run the RANSAC sampling loop
coeff.thDist = 1.5;           % inlier distance threshold; units are in pixels
coeff.source_feat_locs = loc_ft0; % 2D locations of the features in the source t0 image
coeff.target_feat_locs = loc_ft1; % 2D locations of the featuers in the target t1 image
coeff.camera_projective_mat = Pleft;   % intrinsic/projective camera matrix
coeff.randomseed = 2.8;                % use rng(2.8) to set your random seed!!!!!!!
% find those inliers!
% [affT,inlierIdx] = ransac3D(Wt0, Wt1, coeff);
% find those inliers!
[R,t] = ransac3D(Wt0, Wt1, coeff);
%
% EOF
end
function [R_tot,t_tot] = ransac3D(Wt0, Wt1, coeff)
R_tot = zeros(3,3);
t_tot = zeros(3,1);
rng(coeff.randomseed);
maxinlier = 0;
for i = 1:1:coeff.iterNum
%     rsample = randi(min(length(Wt0),length(Wt1)),[1,coeff.numPtsToSample]);
    rsample = randperm(184,3);
%     test = Wt0(rsample(1),:);
    Wt0_samp = [Wt0(rsample(1),:)' Wt0(rsample(2),:)' Wt0(rsample(3),:)'];
    mu_0 = [mean(Wt0_samp(1,:)); mean(Wt0_samp(2,:)); mean(Wt0_samp(3,:))];

    Wt1_samp = [Wt1(rsample(1),:)' Wt1(rsample(2),:)' Wt1(rsample(3),:)'];
    mu_1 = [mean(Wt1_samp(1,:)); mean(Wt1_samp(2,:)); mean(Wt1_samp(3,:))];

    p0_bar = Wt0_samp - mu_0;
    p1_bar = Wt1_samp - mu_1;

    [U, ~, V] = svd(p0_bar * p1_bar');

    R = (V*U')';
    t = mu_0-R*mu_1; %changed this
%     mu_0 = [mean(Wt0_samp(1,:)); mean(Wt0_samp(2,:)); mean(Wt0_samp(3,:)); 1];
%     mu_1 = [mean(Wt1_samp(1,:)); mean(Wt1_samp(2,:)); mean(Wt1_samp(3,:)); 1];
    append_1 = ones([1,length(Wt0)]);
    mu_0_cloud = [Wt0(:,1)'; Wt0(:,2)'; Wt0(:,3)'; append_1];
    mu_1_cloud = [Wt1(:,1)'; Wt1(:,2)'; Wt1(:,3)'; append_1];
    T_trans = [R t; 0 0 0 1];

    rp_xy0 = coeff.camera_projective_mat*inv(T_trans)*mu_0_cloud;
    rp_xy1 = coeff.camera_projective_mat*T_trans*mu_1_cloud;
    
%     rp_xy0_hmg = [rp_xy0(1)/rp_xy0(3); rp_xy0(2)/rp_xy0(3)];
%     rp_xy1_hmg = [rp_xy1(1)/rp_xy1(3); rp_xy1(2)/rp_xy1(3)];
    rp_xy0_hmg = [rp_xy0(1,:)./rp_xy0(3,:); rp_xy0(2,:)./rp_xy0(3,:)];
    rp_xy1_hmg = [rp_xy1(1,:)./rp_xy1(3,:); rp_xy1(2,:)./rp_xy1(3,:)];
    
    Jt0 = coeff.source_feat_locs';
    Jt1 = coeff.target_feat_locs';
    dist_0 = (Jt0(1,:)-rp_xy1_hmg(1,:)).^2+(Jt0(2,:)-rp_xy1_hmg(2,:)).^2;
    dist_1 = (Jt1(1,:)-rp_xy0_hmg(1,:)).^2+(Jt1(2,:)-rp_xy0_hmg(2,:)).^2;
    error = dist_0 + dist_1;
    
    numinlier = 0;
%     for j = 1:1:length(dist_0)
%         if dist_0(j) < coeff.thDist || dist_1(j) < coeff.thDist
%             numinlier = numinlier + 1;
%         end
%     end
    for j = 1:1:length(dist_0)
        if error(j) < coeff.thDist
            numinlier = numinlier + 1;
        end
    end
    
    if numinlier > maxinlier
        maxinlier = numinlier;
        R_tot = R;
        t_tot = t;
    end
end
end
%
function [point_cloud, XL, matchedPointsL_f] = gen_pointcloudandfeatures(Ileft, Iright, Pleft, Pright)
%
% performs image triangulation on stereo pair given by Ileft and Iright

%
% extract and match SURF features from the left and right images
[matchedPointsL, matchedPointsR,matchedPointsL_f,matchedPointsR_f] = detectFeatureMatches(Ileft, Iright); 
XL = matchedPointsL.Location;
XR = matchedPointsR.Location;
% XL = Pleft*XL;
% XR = Pright*XR;
x_L = XL(:,1);
y_L = XL(:,2);
x_R = XR(:,1);
y_R = XR(:,2);

%For all point correspondences
for i = 1:size(matchedPointsL,1)
    % For all of the matched/corresponding points between ImageL and ImageR
    %
    % Perform linear triangulation
    % YOU NEED TO IMPLEMENT THIS FUNCTION
    % point_cloud = linear_triangulation(point_left, point_right);
    %
    %
    A = [x_L(i)*Pleft(3,:)-Pleft(1,:);
         y_L(i)*Pleft(3,:)-Pleft(2,:);
         x_R(i)*Pright(3,:)-Pright(1,:);
         y_R(i)*Pright(3,:)-Pright(2,:)];
    [U,~,V] = svd(A);
%     E = eigs(double(V),4,'smallestabs');
    E = double(V(1:3,4))/double(V(4,4));
    point_cloud(i,:) = E;
%     point_cloud = point_cloud';
%     [U, ~, V] = svd(p1_bar * p2_bar');
%     R = V*U';
%     t = -R*mu1+mu2;
end

%
% Visualization Code
% scatter3(point_cloud(:,1), point_cloud(:,2), point_cloud(:,3))

% EOF
end