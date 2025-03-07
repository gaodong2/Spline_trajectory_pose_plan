function [X, Jp, Ap, Vp, Sp] = P2PMultiAxisDoubleSTrajectory(sampleTime, input)

[rows, columns] = size(input);  %获取参数表的行数与列数
if columns ~= 5
    error('wrong parameters!')
end

T = zeros(1,rows);
N = zeros(1,rows);
Sr = zeros(1,rows);
Vm = zeros(1,rows);
Am = zeros(1,rows);
Jm = zeros(1,rows);
Vj = zeros(1,rows);
Dr = zeros(1,rows);
Da = zeros(1,rows);
Dir = zeros(1,rows);
Dup = zeros(1,rows);
T1 = zeros(1,rows); T2 = zeros(1,rows); T3 = zeros(1,rows); T4 = zeros(1,rows);
T5 = zeros(1,rows); T6 = zeros(1,rows); T7 = zeros(1,rows); Tj = zeros(1,rows);

for i = 1:rows
    Sr(i) = input(i,2) - input(i,1);
    Dr(i) = abs(Sr(i));
    Dir(i) = sign(Sr(i));
    Vm(i) = input(i,3);
    Am(i) = input(i,4);
    Jm(i) = input(i,5);
    Vj(i) = power(Am(i), 2) / (2 * Jm(i));
    if (Vm(i) - 2 * Vj(i) >= 0)
        Tj(i) = Am(i) / Jm(i);
        Da(i) = 2 * power(Am(i), 3) / power(Jm(i), 2);
        Dup(i) = power(Vm(i), 2) / (2 * Am(i)) + Am(i) * Vm(i) / (2 * Jm(i));
        if (Dr(i) >= 2 * Dup(i))
            T1(i) = Tj(i);  T3(i) = Tj(i);  T5(i) = Tj(i);  T7(i) = Tj(i);
            T2(i) = (Vm(i) - 2 * Vj(i)) / Am(i);    T6(i) = T2(i);
            T4(i) = (Dr(i) - 2 * Dup(i)) / Vm(i);
        else
            if (Da(i) <= Dr(i))
                T1(i) = Tj(i);  T3(i) = Tj(i);  T5(i) = Tj(i);  T7(i) = Tj(i);
                T2(i) = (sqrt(power(Am(i), 4) + 4 * Dr(i) * Am(i) * power(Jm(i), 2)) - 3 * power(Am(i), 2)) / (2 *Jm(i) * Am(i));
                T6(i) = T2(i);
                T4(i) = 0;
            else
                T1(i) = power(Dr(i) / (2 * Jm(i)), 1/3);  T3(i) = T1(i);  T5(i) = T1(i);  T7(i) = T1(i);
                T2(i) = 0; T4(i) = 0; T6(i) = 0; 
            end
        end

    else
        Am(i) = sqrt(Vm(i) * Jm(i));
        Tj(i) = Am(i) / Jm(i);
        Da(i) = 2 * power(Am(i), 3) / power(Jm(i), 2);
        Vj(i) = power(Am(i), 2) / (2 * Jm(i));
        T2(i) = 0; T6(i) = 0;
        if (Da(i) <= Dr(i))
            T1(i) = Tj(i);  T3(i) = Tj(i);  T5(i) = Tj(i);  T7(i) = Tj(i);
            T4(i) = (Dr(i) - Da(i)) / Vm(i);
        else
            T1(i) = power(Dr(i) / (2 * Jm(i)), 1/3);  T3(i) = T1(i);  T5(i) = T1(i);  T7(i) = T1(i);
            T4(i) = 0;
        end
    end
    T(i) = 4 * T1(i) + 2* T2(i) + T4(i);
    N(i) = round(T(i) / sampleTime);
    X = 0:sampleTime:(N)*sampleTime;
    Jp = zeros(1,N(i));
    Ap = zeros(1,N(i));
    Vp = zeros(1,N(i));
    Sp = zeros(1,N(i));

    for j = 1: N(i)
        if (j * sampleTime < T1(i))
            Jp(j) = Dir(i) * Jm(i);
        elseif (j * sampleTime >= T1(i) && j * sampleTime < T1(i) + T2(i))
            Jp(j) = 0;
        elseif (j * sampleTime >= T1(i) + T2(i) && j * sampleTime < T1(i) + T2(i) + T3(i))
            Jp(j) = -1 * Dir(i) * Jm(i);
        elseif (j * sampleTime >= T1(i) + T2(i) + T3(i) && j * sampleTime < T1(i) + T2(i) + T3(i) + T4(i))
            Jp(j) = 0;
        elseif (j * sampleTime >= T1(i) + T2(i) + T3(i) + T4(i) && j * sampleTime < T1(i) + T2(i) + T3(i) + T4(i) + T5(i))
            Jp(j) = -1 * Dir(i) * Jm(i);
        elseif (j * sampleTime >= T1(i) + T2(i) + T3(i) + T4(i) + T5(i) && j * sampleTime < T1(i) + T2(i) + T3(i) + T4(i) + T5(i) + T6(i))
            Jp(j) = 0;
        elseif (j * sampleTime >= T1(i) + T2(i) + T3(i) + T4(i) + T6(i))
            Jp(j) = Dir(i) * Jm(i);
        end
    end
    pp = 0;
    pm = 0;
    pmj = [];
    mp = 0;
    mpj = [];
    mm = 0;
    mmj = [];
    for j = 1:N(i)
        if j < N(i)/2
            if Jp(j) > 0
                pp = pp + 1;
            elseif Jp(j) < 0
                pm = pm + 1;
                pmj = [pmj, j];
            end
        else
            if Jp(j) < 0
                mp = mp + 1;
                mpj = [mpj, j];
            elseif Jp(j) > 0
                mm = mm + 1;
                mmj = [mmj, j];
            end   
        end
    end
    
    if(pp~=pm)
        Jp(pmj(end)) = 0; 
    end
    if(mp~=mm)
        Jp(mpj(1)) = 0;
        Jp(mmj(1)) = Jp(mmj(1)-1);
    end
    for j = 1: N(i)
        if j == 1
            Ap(j) = Jp(j) * sampleTime;
            Vp(j) = 1 / 2 * Jp(j) * power(sampleTime, 2);
            Sp(j) = 1 / 6 * Jp(j) * power(sampleTime, 3);
        else
            Ap(j) = Ap(j - 1) + Jp(j) * sampleTime;
            Vp(j) = Vp(j - 1) + Ap(j - 1) * sampleTime + 1 / 2 * Jp(j) * power(sampleTime, 2);
            Sp(j) = Sp(j - 1) + Vp(j - 1) * sampleTime + 1 / 2 * Ap(j - 1) * power(sampleTime, 2) + 1 / 6 * Jp(j) * power(sampleTime, 3);
        end
    end

    Jp = [Jp,0];
    Ap = [Ap,0];
    Vp = [Vp,0];
    Sp = [Sp,Sr(i)];

    figure
    subplot(2,2,1);
    plot(X,Sp);
    title('Position');

    subplot(2,2,2);
    plot(X,Vp);
    title('Speed');

    subplot(2,2,3);
    plot(X,Ap);
    title('Acceleration');

    subplot(2,2,4);
    plot(X,Jp);
    title('Jerk');

end
end