function knots = knots_uniform(n, p)
knots = zeros(1, n + p + 2);
for i = 1:(n + p + 2)
    if i <= p + 1
        knots(i) = 0;
    elseif i > n + 1
        knots(i) = 1;
    else
        knots(i) = (i - p - 1) / (n - p + 1);
    end
end
end
