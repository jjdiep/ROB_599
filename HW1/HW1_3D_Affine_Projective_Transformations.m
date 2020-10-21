clear; close all; clc;
pnt_2D = transform3Dto2D(1000);

function pnt_2D = transform3Dto2D(ws)
pnt = gen_hourglass(100, 10);

theta1 = deg2rad(105);
theta2 = deg2rad(45);
s = 20;
T1 = trans([250; 250; 0]);

T2_rot = rot([0; 0; 0], [0; 1; 0], theta1);
T2_scale = scale([0; 0; 0], eye(3), [50;1;0.25]);
% T2_scale = [50 0 0 0;
%              0 1 0 0;
%              0 0 .25 0;
%              0 0 0 1];
T2 = T2_scale * T2_rot;
% T2 = T2_rot;

t3 = [0;350;0];
T3 = trans(t3(1:3));

T3_A = T3*T2*T1;
% revT3 = [inv(T3_A(1:3,1:3)), -inv(T3_A(1:3,1:3)) * T3_A(1:3,4);...
%         0, 0 , 0, 1];
ref4 = T3_A * [0; 0; 0; 1]; % origin of the current frame resolved in the world frame
% n4_vect = inv(T3*T2_rot*T1)*[1;0;0;0]; % x-axis of the current frame resolved in the world frame
n4_vect = T3*T2_rot*T1*[1;0;0;0]; % x-axis of the current frame resolved in the world frame
n4 = n4_vect/norm(n4_vect);
T4 = rot(ref4(1:3), n4(1:3), theta2);

T4_x = [1 0 0 0;
      0 cos(theta2) -sin(theta2) 0;
      0 sin(theta2) cos(theta2)  0;
      0 0 0 1];
T4_A = T4*T3*T2*T1;
pnt_3D = T4_A * pnt';

p_proj = [1 0 0 0;
          0 1 0 0;
          0 0 1 0;
          0 0 1 0];
       
% T4 = revT3* T4 * T3_A;
% T4_A = T3*T2*T1;
% revT4 = [inv(T4_A(1:3,1:3)), -inv(T4_A(1:3,1:3)) * T4_A(1:3,4);...
%         0, 0 , 0, 1];
% revT4 = inv(T4_A);
% T4 = revT4;
% pnt_3Dv2 = trans([-3235; 600; -60])*T2_rot*T4_x*inv(T2_rot)*trans([3235; -600; 60])*T4_A*pnt';
% T4 = T4;
% T4 = T3_A*T4;

% T = T4 * T3 * T2 * T1;
% T = T4* T3 * T2 * T1;
% pnt_3Dv2 = ceil(T3 * T2 * T1 * pnt');

% pnt_3Dv2 = pnt * T1;
% figure(2)
% clf()
figure(1)
clf()
% scatter3(pnt_3Dv2(1,:),pnt_3Dv2(2,:),pnt_3Dv2(3,:),100,'.m')
hold on;
P = [1 0 0 0;
     0 1 0 0;
     0 0 1 0];
pnt_proj = P * pnt_3D;
for k = 1:size(pnt_proj, 2)
    pnt_proj(:, k) = pnt_proj(:, k) / pnt_proj(3, k);
end
ptrans = [1 0 abs(min(pnt_proj(1,:)));
          0 1 abs(min(pnt_proj(2,:)));
          0 0 1];
pnt_proj = ptrans*pnt_proj;
scaleX = 999/(max(pnt_proj(1,:))-min(pnt_proj(1,:)));
scaleY = 999/(max(pnt_proj(2,:))-min(pnt_proj(2,:)));
pscale = [scaleX 0 0;
          0 scaleY 0;
          0 0 1];
pnt_proj = pscale*pnt_proj;
pnt_2D = pnt_proj;
% for x = 1:size(pnt_proj, 2)
%     pnt_proj(1,x) = pnt_proj(1,x) + abs(min(pnt_proj(1,:)));
% abs(min(pnt_proj(2,:)));
% end

figure(3)
clf
scatter(pnt_proj(1,:),pnt_proj(2,:),'filled','k')
%% Scale and translate `pnt_proj` such that min(x) = min(y) = 0, max(x) = max(y) = 999
% pnt_2D = ???;

figure(1)
% scatter(pnt_2D(1, :), pnt_2D(2, :), 100, '.')
scatter3(pnt_3D(1,1:end-1),pnt_3D(2,1:end-1),pnt_3D(3,1:end-1),s,'filled','b')
scatter3(pnt_3D(1,end),pnt_3D(2,end),pnt_3D(3,end),s,'filled','r')
grid on
% axis([1, ws, 1, ws] - 1)
xlabel('x')
ylabel('y')
zlabel('z')
end
function T = scale(ref, basis, s)
% ref:
%   (3, 1) double
%   reference point
% basis
%   (3, 3) double
%   basis vectors of the frame
% s:
%   (3, 1) double
%   scaling factor along basis vectors
S = [basis * diag(s) * basis', zeros(3, 1); zeros(1, 3), 1];
T = trans(ref) * S * trans(-ref);
end
function T = trans(t)
% t:
%   (3, 1) double
%   displacement vector
T = [eye(3), t; zeros(1, 3), 1];
end
function T = rot(ref, n, theta)
% ref:
%   (3, 1) double
%   a reference point on the rotation axis
% n:
%   (3, 1) double
%   direction of the rotation axis, length must be 1.0
% theta:
%   (1, 1) double
%   rotation angle

% https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
K = [0, -n(3), n(2); n(3), 0, -n(1); -n(2), n(1), 0];
R = eye(3) + sin(theta) * K + (1 - cos(theta)) * K^2;

T = trans(ref) * [R, zeros(3, 1); zeros(1, 3), 1] * trans(-ref);
end