% 定义 B 样条基函数及其导数
function [N0, N1, N2, N3] = bspline_basis(p, knots, u)
    % 计算基函数值
    N0 = zeros(1, length(knots) - p - 1);
    N1 = zeros(1, length(knots) - p - 1);
    N2 = zeros(1, length(knots) - p - 1);
    N3 = zeros(1, length(knots) - p - 1);

    for i = 1:length(N0)
        N0(i) = bspline_basis_value(p, knots, i, u);
        N1(i) = bspline_basis_derivative(p, knots, i, u, 1);
        N2(i) = bspline_basis_derivative(p, knots, i, u, 2);
        N3(i) = bspline_basis_derivative(p, knots, i, u, 3);
    end
end
