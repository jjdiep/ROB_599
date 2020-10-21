velodyne = load('velodyne.mat');
p = velodyne.data;
[R, t] = icp_3d(p{1}, p{2}); 

function [R_total, t_total] = icp_3d(p_t,p_s)
%%
figure(1)
clf()

hold on

%% Get points
N = 500;
p_t = downsample(p_t, N);
p_s = downsample(p_s, N);

scatter3(p_t(1, :), p_t(2, :), p_t(3, :),'bx');
xlabel('x')
ylabel('y')
zlabel('z')

%%
k = 0;
R_total = eye(3);
t_total = zeros(3, 1);
while 1
    %% Update p2, the current source cloud
    p2 = R_total' * (p_s - t_total);
    
    %% Get p1, the matching points from the target cloud
    [match, min_dist] = match_pnt(p_t, p2);
    d = min_dist/std(min_dist);
    w = exp(-d)/sum(exp(-d));
    p1 = p_t(:, match);
    
    %% Plot p2
    if k == 0
        l2 = scatter3(p2(1, :), p2(2, :), p2(3,:),'r+');
        pause(1)
    else
        l2 = scatter3(p2(1, :), p2(2, :), p2(3,:),'r+');
    end
    
    
    %% Centroids of p1 and p2
%     mu1=[mean(p1(1,:));mean(p1(2,:))];
    mu1 = p1 * w;
%     mu2=[mean(p2(1,:));mean(p2(2,:))];
    mu2 = p2 * w;

    %% Center the two clouds
    p1_bar = (p1 - mu1).*w';
    p2_bar = (p2 - mu2).*w';
    
    %% Estimate the rotation and translation
    [U, ~, V] = svd(p1_bar * p2_bar');
    R = (V*U')';
    t = mu1-R*mu2;
    
    %% Update R_total and t_total
    t_total = t_total-(R'*R_total)*t;
    R_total = (R*R_total')';

    %% Terminate when [R, t] is very close to [I, 0]
    delta = norm([R - eye(3), t]);
    if delta > sqrt(eps)
        pause(0.02)
        delete(l2)
        k = k + 1;
        title(k)
    else
        break
    end
end

end

