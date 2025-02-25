function s = sCurve(t, k)
    % 时间缩放的 S 形曲线
    % 输入：
    % t: 插值参数，范围 [0, 1]
    % k: 缩放因子（k > 1 会减小最大角速度）
    % 输出：
    % s: 插值参数，范围 [0, 1]

    t_scaled = t / k; % 时间缩放
    s = 3 * t_scaled.^2 - 2 * t_scaled.^3;
end

