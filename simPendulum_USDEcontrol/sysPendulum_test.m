% Matlab code for simulating a simple pendulum
clc
clear
close all

P.g = 9.81;     %(m/s^2) gravity acceleration
P.l = 1.0;      %(m) length of pendulum
P.c = 1.0;     %(1/s) viscous damping constant
P.m = 1.0;     %(kg) pendulum mass

th0 = 0.5;  %(rad) initial angle of pendulum, from -j axis
w0 = 0.0;   %(rad/s) initial angular rate of pendulum
z0 = [th0;w0];   %Initial state vector

tSpan = [0,20]; %(s) [start, end] times for simulation
dtMax = 0.01;  %(s) maximum allowable time step

%% des
des.th = 0.1; 
des.w = 0.0;
des.a = 0.0;

des.Kp = 100;
des.Kd = 2*sqrt(des.Kp);

des.m_hat = P.m* 0.788 ;% 质量估计
des.g_hat = P.g* 0.888 ;% 重力估计

%%  sim

nStep = ceil(diff(tSpan)/dtMax);  %Number of simulation steps
t = linspace(tSpan(1),tSpan(2),nStep);  %(s) time vector
z = zeros(2,nStep); %  Initialize the state matrix
z(:,1) = z0;
z1 = zeros(2,nStep); %  Initialize the state matrix
z1(:,1) = z0;

tau_prin = zeros(1,nStep); %  用于记录打印的数据
tauf_prin = zeros(1,nStep); %  
P_prin = zeros(1,nStep); %  
Pf_prin = zeros(1,nStep); %  


stopStep1 = ceil(nStep /3); % 第一次改变质量的时间点
stopStep2 = ceil(nStep / 3 * 2); % 第二次

m1 = 10;% 第一次质量变化
m2 = 20;

%Run the simulation, using Euler integration
for i=2:nStep
    dt = t(i)-t(i-1);
    zNow = z(:,i-1);

    if( i > stopStep1)
        P.m = m1;
    end

    if( i > stopStep2)
        P.m = m2;
    end

    % trac desired 期望的sin曲线
    fk = 2*pi*2; % 角频率

    traj.th(i) = 0.1 * sin(fk*i*dtMax); % 期望角度
    traj.w(i) = fk*0.1 * cos(fk*i*dtMax); % 期望角速度
    traj.a(i) = -fk*fk*0.1 * sin(fk*i*dtMax); % 期望角加速度


    des.th = traj.th(i);
    des.w = traj.w(i);
    des.a = traj.a(i);

    tau = control(zNow, des); % PD+G 控制
    z(:,i) = zNow + dt*dynamics(zNow,P, tau);% i-1
end

%Run the simulation, using USDE control
for i=2:nStep
    dt = t(i)-t(i-1);
    zNow = z1(:,i-1);

    if( i > stopStep1)
        P.m = m1;
    end

    if( i > stopStep2)
        P.m = m2;
    end

    % trac
    fk = 2*pi*2;

    traj.th(i) = 0.1 * sin(fk*i*dtMax);
    traj.w(i) = fk*0.1 * cos(fk*i*dtMax);
    traj.a(i) = -fk*fk*0.1 * sin(fk*i*dtMax);


    des.th = traj.th(i);
    des.w = traj.w(i);
    des.a = traj.a(i);

    if i==2
        % reset
        last.P = 0;
        last.tau = 0;
        last.tau_last = 0;
        last.g_last = 0;
    end

    [tau, now, print] = USDEcontrol(zNow, des, last); % USDE控制


    % memeory
    last.P = now.P; % 记录控制量用于下次计算
    last.tau = tau;
    last.tau_last = last.tau;
    last.g_last = now.g_hat;

    % save
    tau_prin(i) = tau;
    tauf_prin(i) = print.tauf;
    P_prin(i) = now.P;
    Pf_prin(i) = print.Pf;

    z1(:,i) = zNow + dt*dynamics(zNow,P, tau);% i-1
end

%Generate a plot of the result:
figure(1); clf;
subplot(311); 
plot(t,traj.th,'r');
hold on
plot(t,z(1,:),'g');
hold on
plot(t,z1(1,:),'b');
xlabel('Time (s)'); ylabel('Angle (rad)');
title('Simple Damped Pendulum');
grid on 
legend('desired','PDcontrol','USDEcontrol');



subplot(312); 
plot(t,traj.w,'r');
hold on
plot(t,z(2,:),'g');
hold on
plot(t,z1(2,:),'b');
xlabel('Time (s)'); ylabel('Rate (rad/s)');
grid on 
subplot(313); 
plot(t,z(1,:)-traj.th(1,:));
hold on
plot(t,z1(1,:)-traj.th(1,:));
xlabel('Time (s)'); ylabel('e (rad)');
grid on 

% figure
% plot(t,P_prin,'r');
% hold on
% plot(t,Pf_prin,'b');


% end

function dz = dynamics(z,P, tau)
% Compute the dynamics of a simple damped pendulum:
th = z(1,:); w = z(2,:); %Unpack the state vector
dth = w;
dw = -(P.g/P.l)*sin(th) - (P.c/(P.m*P.l*P.l))*w + tau /(P.m*P.l*P.l) ;% 随机噪声+ 0.1*randn(size(th))
dz = [dth; dw];  %Pack up state derivative
end


function tau = control(z, des)
% control it
th = z(1,:); w = z(2,:); % current state
dth = w;
g_hat = des.m_hat * des.g_hat * sin(th);
tau = 1 / des.m_hat * (des.Kp * (des.th - th) + des.Kd * (des.w - dth) + g_hat) ;
end

% USDE Control
function [tau,now,print] = USDEcontrol(z, des, last)
% control it
th = z(1,:); w = z(2,:); % current state
dth = w;

% est
P = des.m_hat * dth; % M*qdot 动量

% filter
Pf = online_Lowpassfilter(last.P, 0.01, 5, P);
tauf = online_Lowpassfilter(last.tau_last, 0.01, 5, last.tau);

k = 0.001; % 一阶观测器系数
d_hat = (P - Pf) / k + des.m_hat * des.g_hat * sin(th) - tauf; % 估计的扰动

% con
e = (des.th - th);
edot = (des.w - dth);

n = 5;
S = edot + n * e;

K = 100;
rho = des.w + n * e;
rhodot = des.a + n * edot;

g_hat = des.m_hat * des.g_hat * sin(th);
g_hatf = online_Lowpassfilter(last.g_last, 0.01, 5, g_hat);

tau = K * S + des.m_hat * rhodot + 0 * rho + g_hatf - d_hat;

print.tauf = tauf;
print.Pf = Pf;

now.P = P;
now.g_hat = g_hat;

end

function output_data= online_Lowpassfilter(last,dt,cutfreq, input_data) %在线一阶低通滤波

output_data = last + dt * (input_data - last) * 2 * pi * cutfreq;

end