function [s, v, t] = trapezoidal_velocity_profile(v0, vt, v_max, a_max, s)
    % 参数设置
    % 初始速度 (m/s)
    % 目标速度 (m/s)
    % 最大加速度 (m/s^2)
    d_max = a_max;       % 最大减速度 (m/s^2)
    % 总位移 (m)
    % 最大速度 (m/s)
    dt = 0.01;       % 时间步长 (s)

% v0 = 0.0206
% vt = 0
% v_max = 0.0206
% a_max = 1
% s = 0.1091


    % 计算加速和减速时间
    t_a = (v_max - v0) / a_max;  % 加速时间
    t_d = (v_max - vt) / d_max;  % 减速时间

    % 计算加速和减速距离
    s_a = v0 * t_a + 0.5 * a_max * t_a^2;  % 加速距离
    s_d = v_max * t_d - 0.5 * d_max * t_d^2;  % 减速距离

    % 计算匀速阶段
    s_c = s - s_a - s_d;  % 匀速距离
    if s_c < 0
        % 如果无法达到最大速度，重新计算
        v_max = sqrt((2 * a_max * d_max * s + d_max * v0^2 + a_max * vt^2) / (a_max + d_max));
        t_a = (v_max - v0) / a_max;
        t_d = (v_max - vt) / d_max;
        s_a = v0 * t_a + 0.5 * a_max * t_a^2;
        s_d = v_max * t_d - 0.5 * d_max * t_d^2;
        s_c = 0;  % 没有匀速阶段
    end
    t_c = s_c / v_max;  % 匀速时间

    % 总时间
    total_time = t_a + t_c + t_d;

    % 生成时间序列
%     t = 0:dt:total_time;
    t = dt:dt:total_time;
    % 初始化速度数组
    v = zeros(size(t));
    s = zeros(size(t));
    % 计算速度曲线
    for i = 1:length(t)
        if t(i) < t_a
            % 加速阶段
            v(i) = v0 + a_max * t(i);
            s(i) = v0 * t(i) + 1/2 * a_max * t(i)^2;
        elseif t(i) < t_a + t_c
            % 匀速阶段
            v(i) = v_max;
            s(i) = v_max * (t(i) - t_a) + s_a;
        else
            % 减速阶段
            v(i) = v_max - d_max * (t(i) - t_a - t_c);
            s(i) = v_max * (t(i) - t_a - t_c) - 1/2 * d_max * (t(i) - t_a - t_c)^2 + s_a  + s_c;
        end



    end
        t = [t, total_time];
        v = [v, vt];
        s = [s, s_a + s_d + s_c];
    % 绘制速度曲线
%     figure;
%     subplot(2, 1, 1);
%     plot(t, v, 'b', 'LineWidth', 2);
%     xlabel('时间 (s)');
%     ylabel('速度 (m/s)');
%     title('梯形速度规划（初始速度和目标速度不为零且不同）');
%     grid on;
% 
%     subplot(2, 1, 2);
%     plot(t, s, 'b', 'LineWidth', 2);
%     xlabel('时间 (s)');
%     ylabel('位移 (m)');
%     title('梯形速度规划（初始速度和目标速度不为零且不同）');
%     grid on;
end