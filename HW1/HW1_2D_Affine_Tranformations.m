% HW1-2D Affine Transformations of Apple Image
close all; clear; clc;
new_apple = apple_transformations('cartoonapple.jpeg');
figure; imshow(new_apple);

function newapple = apple_transformations(apple_filepath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

apple = load_apple(apple_filepath);
apple = uint8(apple);
world = zeros(500,'uint8');
size(world);
% translation
% rotation
theta = deg2rad(-90);
theta2 = deg2rad(45);

%% Second Transformation Algo
T_trans1 = [1 0 250; 
            0 1 200; 
            0 0 1];
T_rot1 = [cos(theta) -sin(theta) 0; 
          sin(theta) cos(theta) 0; 
          0 0 1];
T_trans2 = [1 0 0; 
            0 1 100; 
            0 0 1];
T_rot2 = [cos(theta2) -sin(theta2) 0; 
          sin(theta2) cos(theta2) 0;
          0 0 1];
T_scale1 = [2 0 0;
            0 2 0;
            0 0 1];
T_trans3 = [1 0 0; 
            0 1 -150; 
            0 0 1];
T_A = T_trans3*T_rot2*T_trans1*T_rot1*T_trans2*T_scale1;

revT = [inv(T_A(1:2,1:2)), -inv(T_A(1:2,1:2)) * T_A(1:2,3);...
        0, 0 ,1];
for i = 0:1:499
    for j = 0:1:499
        p_old = round(revT*[i;j;1]);
        if p_old(1) >= 0 && p_old(2) >= 0
            if p_old(1) < 100 && p_old(2) < 100
                world(j+1,i+1) = apple(p_old(2)+1,p_old(1)+1);
    %           p_new = ceil(T_trans3*T_trans2*T_trans1*T_rot2*T_rot1*T_scale1*p);
            end
        else
            world(j+1,i+1) = 0;
        end
    end
end

newapple = world;


end

