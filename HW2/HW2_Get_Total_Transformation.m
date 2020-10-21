[R_total, t_total] = get_total_transformation(R, t)

function [R_total, t_total] = get_total_transformation(R, t)
    R_total = cell(length(R), 1);
    t_total = cell(length(R), 1);
    C_total = cell(length(R), 1);
    
    R_total{1} = R{1};
    t_total{1} = t{1};
    C_total{1} = [R_total{1} t_total{1}; 0 0 0 1];
    for i = 2:1:length(R)
          Curr = [R{i} t{i}; 0 0 0 1];
          C_total{i} = Curr * C_total{i-1};
          R_total{i} = C_total{i}(1:3,1:3);
          t_total{i} = C_total{i}(1:3,4);
    end
end



