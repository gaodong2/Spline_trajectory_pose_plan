function q_inv = quatInverse(q)
    % 四元数的逆
    q_inv = [q(1), -q(2), -q(3), -q(4)];
end

