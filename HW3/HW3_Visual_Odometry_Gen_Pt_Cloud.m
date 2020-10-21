% load in images
Ileft = rgb2gray(imread('It1_L.png'))
Iright = rgb2gray(imread('It1_R.png'))
% load in calibration data
calibname = 'calib.txt';
T = readtable(calibname, 'Delimiter', 'space', 'ReadRowNames', true, 'ReadVariableNames', false);
B = table2array(T);
Pleft = vertcat(B(1,1:4), B(1,5:8), B(1,9:12)); % left camera [3x4] projection matrix
Pright = vertcat(B(2,1:4), B(2,5:8), B(2,9:12)); % right camera [3x4] projection matrix
% call function!
[point_cloud] = gen_pointcloud(Ileft, Iright, Pleft, Pright)
%

function point_cloud = gen_pointcloud(Ileft, Iright, Pleft, Pright)
%
% performs image triangulation on stereo pair given by Ileft and Iright

%
% extract and match SURF features from the left and right images
[matchedPointsL, matchedPointsR] = detectFeatureMatches(Ileft, Iright); 
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
scatter3(point_cloud(:,1), point_cloud(:,2), point_cloud(:,3))

% EOF
end