function s = sCurveQuintic(t)
    % 五次多项式 S 形曲线
    % 输入：
    % t: 插值参数，范围 [0, 1]
    % 输出：
    % s: 插值参数，范围 [0, 1]

    s = 10 * t.^3 - 15 * t.^4 + 6 * t.^5;
end

