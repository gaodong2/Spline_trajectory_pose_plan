function  q_out = do_slerp(q1,q2,t)
% slerp插值
EPS = 1e-9;
C = dot(q1,q2);

if ((1 - C) <= EPS) 
    % 两姿态点积接近1，即夹角非常接近于０，采用直线插补
    q_out=q1*(1-t)+q2*t; % avoiding divisions by number close to 0
    q_out = quatnormalize(q_out);
    return;
elseif((1 + C) <= EPS) 
    % 两姿态点积接近-1，即夹角非常接近于18０
    qtemp(1) = q2(4); qtemp(2) = -q2(3); qtemp(3)= q2(2); qtemp(4) = -q2(1); % rotating one of the unit quaternions by 90 degrees -> q2
    q2 =qtemp;
end
% 插补
q1_inv    = quatinv(q1);
q1_inv_q2 = quatmultiply(q1_inv,q2);
omega     = quatLog(q1_inv_q2);
q_out = quatmultiply(q1,quatexp(omega*t));
%q_out = quatnormalize(q_out);
end
