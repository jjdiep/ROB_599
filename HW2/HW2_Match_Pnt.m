[p_t, p_s] = get_pnt_2d();
[match, min_dist] = match_pnt(p_t, p_s);

function[match,min_dist] = match_pnt(p_t,p_s)
    for i=1:length(p_s)
        for j=1:length(p_t)
            normdiff(j) = norm(p_t(:,j)-p_s(:,i));
        end
        [minval,ind] = min(normdiff);
        match(i) = ind;
        min_dist(i) = minval;
    end
    match = match';
    min_dist = min_dist';
end