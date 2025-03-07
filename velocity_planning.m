function v = velocity_planning(s_start, s_end, v_max, a_sup, v0, vt)
    % 输入参数
    s = [s_start; s_end];  % 曲线上的位置 (m)
%     v_max = v_sup * ones(size(s));  % 各点的最大速度 (m/s)
    a_max = a_sup * ones(size(s));  % 各点的最大加速度 (m/s^2)
%     v0 = 0;  % 初始速度 (m/s)
%     vt = 0;  % 终止速度 (m/s)

    % 前向传播
    v_forward = zeros(size(s));
    v_forward(1) = v0;
    for i = 2:length(s)
        delta_s = s(i) - s(i-1);
        v_forward(i) = min(v_max(i), sqrt(v_forward(i-1)^2 + 2 * a_max(i) * delta_s));
    end
    
    % 反向传播
    v_backward = zeros(size(s));
    v_backward(end) = vt;
    for i = length(s)-1:-1:1
        delta_s = s(i+1) - s(i);
        v_backward(i) = min(v_max(i), sqrt(v_backward(i+1)^2 + 2 * a_max(i) * delta_s));
    end
    % 可达速度
    v = min(v_forward, v_backward);

%      绘制结果
    figure;
    plot(s, v_max, 'r--', 'LineWidth', 2); hold on;
    plot(s, v_forward, 'g-.', 'LineWidth', 2);
    plot(s, v_backward, 'b-.', 'LineWidth', 2);
    plot(s, v, 'k', 'LineWidth', 2);
    xlabel('位置 s (m)');
    ylabel('速度 v (m/s)');
    legend('最大速度限制', '前向传播速度', '反向传播速度', '可达速度');
    title('速度规划');
    grid on;
end
