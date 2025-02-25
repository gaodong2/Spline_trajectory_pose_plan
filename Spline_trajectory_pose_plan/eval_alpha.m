function alpha = eval_alpha(s,i,L)
% 全局细分到局部细分
k = s*(L-1)+1;

if((i>=k)&&(i<k+1))
    alpha=k-(i-1);
else
    alpha=0;
end

end

