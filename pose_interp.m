  pose = [ 0.001/180*pi, -0.001/180*pi, 90.000/180*pi
-60.391/180*pi, -0.920/180*pi, 88.689/180*pi;
-60.401/180*pi,  0.149/180*pi, 90.215/180*pi;
 0.528/180*pi, -0.962/180*pi, 87.596/180*pi];

  q= [];
for i = 1:4
qq = eul2quat(pose(i,:), "ZYX");
q = [q; qq];
end
num_samples = 1000;
t = linspace(0, 1, num_samples);
s = sCurve(t);
VAL = [q(1,:)];
for i =2:999
    val = quat_squad(q',s(i))
    VAL = [VAL; val];
end
VAL = [VAL; q(4,:)];
dt = ones(1000,1);
dt = dt * 3.80038018406146 / 1000;
% dt(1) = 0.05;
% dt(1000) = 0.05;
T = 0;
[omega, alpha] = computeAngularVelocityAndAcceleration(VAL, dt);
for i = 1:999
    T = [T,sum(dt(1:i))];
end
X = 0.001:0.001:0.001*1000;
plot(T, omega(:,1));
hold on;
plot(T, omega(:,2));
hold on;
plot(T, omega(:,3));













