% 计算 B 样条基函数值
function val = bspline_basis_value(p, knots, i, u)
    if p == 0
        val = (knots(i) <= u) && (u < knots(i + 1));
    else
        denom1 = knots(i + p) - knots(i);
        denom2 = knots(i + p + 1) - knots(i + 1);
        
        term1 = 0;
        term2 = 0;
        
        if denom1 ~= 0
            term1 = (u - knots(i)) / denom1 * bspline_basis_value(p - 1, knots, i, u);
        end
        
        if denom2 ~= 0
            term2 = (knots(i + p + 1) - u) / denom2 * bspline_basis_value(p - 1, knots, i + 1, u);
        end
        
        val = term1 + term2;
    end
end

