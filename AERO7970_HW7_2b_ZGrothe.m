% Zane Grothe
% AERO 7970
% HW 7
% 10/21/22

clear all
close all
clc

% Problem 2 ~~~~~~~~~~~~~~~~~~~~

% Part b)

%  Initial conditions
t0 = 0;                           % starting time
tf = 2*pi/1i;                     % final time
n = 500;                          % number of steps
dt = tf/n;                        % step size
t = zeros(1,n+1);                 % time range of result (imag numbers)
for j = 1:(length(t)-1)
    t(j+1) = t(j) + dt;
end
%t = t0:dt:tf;                    % time range of result (real numbers)
z = zeros(1,length(t));           % preallocation
z(1) = 1+1i;                      % initial conditions
F_tz1 = @(t,z) -z/abs(z)^3;       % dz1/dt
F_tz2 = @(t,z) z;                 % dz2/dt

% Color matrix to pick from for plotting
C=[1,0,0; 0,1,0; 0,0,1; .929,.694,.125; 0,1,1; 1,0,1; 0,0,0; .85,.325,.098];
% [ red ; green; blue ; gold          ; cyan ; mag. ; black; brown        ]

%% Forward Euler (red)
% y(t_i+1) = y(t_i) + dt*f(t_i, y(t_i))

% Find velocity
for j = 1:(length(t)-1)
    z(j+1) = z(j) + dt*F_tz1(t(j),z(j));
end

% Find position
for j = 1:(length(t)-1)
    z(j+1) = z(j) + dt*F_tz2(t(j),z(j));
end

% Separate x and y components
z_feuler_x = real(z);
z_feuler_y = imag(z);

% Plot
figure(1)
plot(z_feuler_x,z_feuler_y,'color',C(1,:))
axis square
hold on

%% RK2 (green)
% y(t_i+1) = y(t_i) + dt/2*( k1 + k2 )

% Find velocity
for j = 1:(length(t)-1)
    k1 = F_tz1(t(j),z(j));
    k2 = F_tz1(t(j+1),z(j)+dt*k1);

    z(j+1) = z(j) + dt/2*(k1+k2);
end

% Find position
for j = 1:(length(t)-1)
    k1 = F_tz2(t(j),z(j));
    k2 = F_tz2(t(j+1),z(j)+dt*k1);

    z(j+1) = z(j) + dt/2*(k1+k2);
end

% Separate x and y components
z_RK2_x = real(z);
z_RK2_y = imag(z);

% Plot
figure(1)
plot(z_RK2_x,z_RK2_y,'color',C(2,:))
axis square
hold on

%% RK4 (blue)
% y(t_i+1) = y(t_i) + dt/6*( k1 + 2*k2 + 2*k3 + k4)

% Find velocity
for j = 1:(length(t)-1)
    k1 = F_tz1(t(j),z(j));
    k2 = F_tz1(t(j)+dt/2,z(j)+dt/2*k1);
    k3 = F_tz1(t(j)+dt/2,z(j)+dt/2*k2);
    k4 = F_tz1(t(j)+dt,z(j)+dt*k3);
    
    z(j+1) = z(j) + dt/6*(k1+2*k2+2*k3+k4);
end

% Find position
for j = 1:(length(t)-1)
    k1 = F_tz2(t(j),z(j));
    k2 = F_tz2(t(j)+dt/2,z(j)+dt/2*k1);
    k3 = F_tz2(t(j)+dt/2,z(j)+dt/2*k2);
    k4 = F_tz2(t(j)+dt,z(j)+dt*k3);
    
    z(j+1) = z(j) + dt/6*(k1+2*k2+2*k3+k4);
end

% Separate x and y components
z_RK4_x = real(z);
z_RK4_y = imag(z);

% Plot
figure(1)
plot(z_RK4_x,z_RK4_y,'color',C(3,:))
axis square
hold on
xlabel('x')
ylabel('y')
title('Numerical Solutions for a Two Dimensional Solar System')
legend({'FEuler','RK2','RK4'},'Location','best')
