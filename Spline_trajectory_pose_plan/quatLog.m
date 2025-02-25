function ln_q = quatLog(q)
    % 四元数的对数映射
    eps
    theta = acos(q(1));
    if abs(theta) < eps % 如果 theta 接近 0
        ln_q = [0, 0, 0, 0]; % 返回零向量
    else
        ln_q = [0, q(2:4)] * (theta / sin(theta));
    end
end
