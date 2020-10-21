[R, t] = icp_2d_weighted();

function [R_total, t_total] = icp_2d_weighted()
%%
figure(1)
clf()

hold on
axis equal
axis([-2, 2, -2, 1])

%% Get points
[p_t, p_s] = get_pnt_2d();

plot(p_t(1, :), p_t(2, :), 'bx');
xlabel('x')
ylabel('y')


%%
k = 0;
R_total = eye(2);
t_total = zeros(2, 1);
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
        l2 = plot([p1(1, :); p2(1, :)], [p1(2, :); p2(2, :)], 'k', ...
                p2(1, :), p2(2, :), 'r+');
        pause(1)
    else
        l2 = plot(p2(1, :), p2(2, :), 'r+');
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
    delta = norm([R - eye(2), t]);
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
