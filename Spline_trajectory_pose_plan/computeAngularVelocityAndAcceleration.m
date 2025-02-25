function [omega, alpha] = computeAngularVelocityAndAcceleration(q, dt)
    % 输入：
    % q: Nx4 矩阵，N 个时间步的四元数（每行是一个四元数 [w, x, y, z]）
    % dt: 时间步长
    % 输出：
    % omega: Nx3 矩阵，角速度（每行是 [ωx, ωy, ωz]）
    % alpha: Nx3 矩阵，角加速度（每行是 [αx, αy, αz]）

    N = size(q, 1); % 时间步数
    omega = zeros(N, 3); % 初始化角速度
    alpha = zeros(N, 3); % 初始化角加速度

    % 计算四元数导数 dq/dt
    dq = zeros(N, 4);
    dq(1, :) = (q(2, :) - q(1, :)) / (dt(1));
    for i = 2:N-1
        dq(i, :) = (q(i+1, :) - q(i-1, :)) / (2 * dt(i-1)); % 中心差分
    end
    dq(N, :) = (q(N, :) - q(N-1, :)) / (dt(N-1));
    % 计算角速度 omega
    for i = 1:N
        q_inv = quatInverse(q(i, :)); % 四元数的逆
        omega_q = 2 * quatMultiply(q_inv, dq(i, :)); % 角速度（四元数形式）
        omega(i, :) = omega_q(2:4); % 提取虚部（角速度向量）
    end

    % 计算角加速度 alpha
    for i = 2:N-1
        alpha(i, :) = (omega(i+1, :) - omega(i-1, :)) / (2 * dt(i-1)); % 中心差分
    end
end
