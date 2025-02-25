% 计算 B 样条基函数的导数
function val = bspline_basis_derivative(p, knots, i, u, order)
    if order == 0
        val = bspline_basis_value(p, knots, i, u);
    else
        denom1 = knots(i + p) - knots(i);
        denom2 = knots(i + p + 1) - knots(i + 1);
        
        term1 = 0;
        term2 = 0;
        
        if denom1 ~= 0
            term1 = p / denom1 * bspline_basis_derivative(p - 1, knots, i, u, order - 1);
        end
        
        if denom2 ~= 0
            term2 = -p / denom2 * bspline_basis_derivative(p - 1, knots, i + 1, u, order - 1);
        end
        
        val = term1 + term2;
    end
end