% function s_speed_up()
    % 参数设置
    v0 = 0;     % 初始速度 (m/s)
    vf = 0.3;     % 终止速度 (m/s)
    amax = 0.1;   % 最大加速度 (m/s^2)
    jmax = 0.02;   % 最大加加速度 (m/s^3)
    s = 0.3;     % 总位移 (m)
    
    vj = amax^2 / (2 * jmax);

    if vf >= 2 * vj
    % 有匀加速(有1，2，3)
        t1 = amax / jmax; 
        t3 = t1;

        v1 = v0 + 1/2 * jmax * (t1)^2;
        v3 = vf - 1/2 * jmax * (t3)^2;
        
        t2 = 2*deltas / (v1 + v3);
        
        
        % 总时间
        total_time = t1 + t2 + t3;
        
        % 时间向量
        t = linspace(0, total_time, 1000);
        
        % 初始化速度、加速度和位移
        v = zeros(size(t));
        a = zeros(size(t));
        pos = zeros(size(t));

        % 分段计算速度、加速度和位移
        for i = 1:length(t)
            if t(i) < t1
                % 加加速度阶段
                a(i) = jmax * t(i);
                v(i) = v0 + 0.5 * jmax * t(i)^2;
                pos(i) = v0 * t(i) + (1/6) * jmax * t(i)^3;
            elseif t(i) < t1 + t2
                % 匀加速阶段
                a(i) = amax;
                v(i) = v0 + amax * (t(i) - t1) + 0.5 * jmax * t1^2;
                pos(i) = v0 * t(i) + 0.5 * amax * (t(i) - t1)^2 + (1/6) * jmax * t1^3;
            elseif t(i) <= t1 + t2 + t3
                % 减加速度阶段
                a(i) = amax - jmax * (t(i) - t1 - t2);
                v(i) = vmax - 0.5 * jmax * (t(i) - t1 - t2)^2;
                pos(i) = s_acc + vmax * (t(i) - t1 - t2) - (1/6) * jmax * (t(i) - t1 - t2)^3;
            end
        end


    elseif vf < 2 * vj
    % 无匀加速（有1，3）
        amax = sqrt(vf * jmax);
        t1 = amax / jmax;
        v1 = 1/2 * jmax * t1^2;

        t3 = t1;
        s1 = 1/6 * jmax * t1^3;
        s3 = v1 * t3 + 1/2 * amax * t3^2 - (1/6) * jmax * t3^3;
        s4 = s - s1 - s3;
        % 总时间
        t4 = s4 / vf;
        total_time = t1 + t3 + t4;
        
        % 时间向量
        t = linspace(0, total_time, 1000);
        
        % 初始化速度、加速度和位移
        v = zeros(size(t));
        a = zeros(size(t));
        pos = zeros(size(t));

        % 分段计算速度、加速度和位移
        for i = 1:length(t)
            if t(i) < t1
                % 加加速度阶段
                a(i) = jmax * t(i);
                v(i) = v0 + 0.5 * jmax * t(i)^2;
                pos(i) = v0 * t(i) + (1/6) * jmax * t(i)^3;
            elseif t(i) <= t1 + t3
                % 减加速度阶段
                a(i) = amax - jmax * (t(i) - t1);
                v(i) = v1 + amax * (t(i) - t1) - 0.5 * jmax * (t(i) - t1)^2;
                pos(i) = s1 + v1 * (t(i) - t1) + 1/2 * amax * (t(i) - t1)^2 - (1/6) * jmax * (t(i) - t1)^3;
            elseif t(i) <= t1 + t3 + t4
                % 减加速度阶段
                a(i) = 0;
                v(i) = vf;
                pos(i) = s1 + s3 + vf * (t(i) - t1 -t3);            
            end
        end

%     elseif deltas > s2
%     % 有匀速阶段（有1，2，3，4）
%         s1 = v0 * t1 + (1/6) * jmax * t1^3; % 加加速阶段位移
%         s2 = (v0 + vf) /2 * t2;
%         s3 = vf * t3 - (1/6) * jmax * t3^3;
%         s4 = s - s1 - s2 - s3;    
%         
%         % 总时间
%         total_time = t1 + t2 + t3 + t4;
%         
%         % 时间向量
%         t = linspace(0, total_time, 1000);
%         
%         % 初始化速度、加速度和位移
%         v = zeros(size(t));
%         a = zeros(size(t));
%         pos = zeros(size(t));
%         
%         % 分段计算速度、加速度和位移
%         for i = 1:length(t)
%             if t(i) < t1
%                 % 加加速度阶段
%                 a(i) = jmax * t(i);
%                 v(i) = v0 + 0.5 * jmax * t(i)^2;
%                 pos(i) = v0 * t(i) + (1/6) * jmax * t(i)^3;
%             elseif t(i) < t1 + t2
%                 % 匀加速阶段
%                 a(i) = amax;
%                 v(i) = v0 + amax * (t(i) - t1) + 0.5 * jmax * t1^2;
%                 pos(i) = v0 * t(i) + 0.5 * amax * (t(i) - t1)^2 + (1/6) * jmax * t1^3;
%             elseif t(i) < t1 + t2 + t3
%                 % 减加速度阶段
%                 a(i) = amax - jmax * (t(i) - t1 - t2);
%                 v(i) = vmax - 0.5 * jmax * (t(i) - t1 - t2)^2;
%                 pos(i) = s_acc + vmax * (t(i) - t1 - t2) - (1/6) * jmax * (t(i) - t1 - t2)^3;
%             elseif t(i) <= t1 + t2 + t3 + t4
%                 % 减加速度阶段
%                 a(i) = 0;
%                 v(i) = vf;
%                 pos(i) = s1 + s2 + s3 + vf * (t(i) - t1 - t2 - t3);
%             end
%         end
    end

    % 绘图
    figure;
    subplot(3,1,1);
    plot(t, pos, 'b', 'LineWidth', 1.5);
    xlabel('时间 (s)');
    ylabel('位移 (m)');
    title('位移-时间曲线');
    grid on;
    
    subplot(3,1,2);
    plot(t, v, 'r', 'LineWidth', 1.5);
    xlabel('时间 (s)');
    ylabel('速度 (m/s)');
    title('速度-时间曲线');
    grid on;
    
    subplot(3,1,3);
    plot(t, a, 'g', 'LineWidth', 1.5);
    xlabel('时间 (s)');
    ylabel('加速度 (m/s^2)');
    title('加速度-时间曲线');
    grid on;
% end