function q = quatExp(v)
    % 四元数的指数映射
    theta = norm(v);
    if theta < eps
        q = [1, 0, 0, 0];
    else
        q = [cos(theta), v(2:4) / theta * sin(theta)];
    end
end

