function [val] =  quat_squad(q,s)

L = size(q,2);
val = q(:,1)';

% 若点积小于0进行取负，确保两姿态之间取最短路径
for j=2:L
    C = dot(q(:,j-1),q(:,j));
    if(C<0)
        q(:,j) = -q(:,j);
    end 
end

% 若时间在0,1则直接返回起点和终点
if s==0 
    val=q(:,1);
    return;
elseif s==1
    val=q(:,end);
    return;
end

for j =2:L
    % 全局细分映射到局部细分
    alpha=eval_alpha(s,j,L);
    t= alpha;
    
    if(alpha>0)
        EPS = 1e-9;

        % 计算两个姿态的点积，结果范围[-1,1]
        C = dot(q(:,j-1),q(:,j));

        if ((1 - C) <= EPS) 
            % 姿态之间过于接近采用线性插补
            val=q(:,j-1)'*(1-s)+q(:,j)'*s; 
            val = quatnormalize(val);
	    return;
		end

        if((1 + C) <= EPS)
            % 当姿态夹角接近180，无最短路径，结果不确定
            % 将角度旋转90
            qtemp(1) = q(4,j); qtemp(2) = -q(3,j); qtemp(3)= q(2,j); qtemp(4) = -q(1,j);
            q(:,j) = qtemp';
        end

        % 计算中间点
        qa = get_intermediate_control_point(j-1,q);
        qap1 = get_intermediate_control_point(j,q);

        % 插补
        qtemp1 = do_slerp(q(:,j-1)', q(:,j)', t);
        qtemp2 = do_slerp(qa, qap1, t);
        squad = do_slerp(qtemp1, qtemp2, 2*t*(1-t));
        val = squad;val = quatnormalize(val);
	return;
    end

end

end