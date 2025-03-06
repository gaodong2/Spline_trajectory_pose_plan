function [qa] = get_intermediate_control_point(j,q)
% 插补中间点
% 若点为起点和终点，则直接返回当前值
L = size(q,2);
if(j==1)
    qa =  q(:,1)';
    return;
elseif(j==L)
    qa =  q(:,L)';
    return;
else
    qji=quatinv(q(:,j)');
    qiqm1=quatmultiply(qji,q(:,j-1)');
    qiqp1=quatmultiply(qji,q(:,j+1)');

    ang_vel =-((quatlog(qiqp1)+quatlog(qiqm1))/4); 
    
    qa = quatmultiply(q(:,j)',quatExp(ang_vel));
    %qa = quatnormalize(qa);
end
end
