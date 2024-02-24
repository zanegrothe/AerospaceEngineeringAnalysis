% Zane Grothe
% AERO 7970
% HW 7
% 10/21/22

clear all
close all
clc

% Problem 2 ~~~~~~~~~~~~~~~~~~~~

% Part c)

%  Initial conditions
error = zeros(3,4);                   % perallocation (error)
t0 = 0;                               % starting time
tf = 2*pi/1i;                         % final time
nm = [500,1000,2000,4000];            % number of steps (matrix)
k = 1;                                % counter
while k < 5
    n = nm(k);
    
    dt = tf/n;                        % step size
    t = zeros(1,n+1);                 % time range of result (imag numbers)
    for j = 1:(length(t)-1)
        t(j+1) = t(j) + dt;
    end
    %t = t0:dt:tf;                    % time range of result (real numbers)
    z = zeros(1,length(t));           % preallocation (position)
    z(1) = 1+1i;                      % initial conditions
    F_tz1 = @(t,z) -z/abs(z)^3;       % dz1/dt
    F_tz2 = @(t,z) z;                 % dz2/dt

    % Color matrix to pick from for plotting
    C=[1,0,0; 0,1,0; 0,0,1; .929,.694,.125; 0,1,1; 1,0,1; 0,0,0; .85,.325,.098];
    % [ red ; green; blue ; gold          ; cyan ; mag. ; black; brown        ]

    %% Forward Euler (red)
    % y(t_i+1) = y(t_i) + dt*f(t_i, y(t_i))

    for j = 1:(length(t)-1)
        z(j+1) = z(j) + dt*F_tz1(t(j),z(j));
    end

    for j = 1:(length(t)-1)
        z(j+1) = z(j) + dt*F_tz2(t(j),z(j));
    end

    error(1,k) = abs(z(1)-z(end));

    figure(1)
    loglog(error(1,k),n,'o','color',C(k,:))
    hold on

    %% RK2 (green)
    % y(t_i+1) = y(t_i) + dt/2*( k1 + k2 )

    for j = 1:(length(t)-1)
        k1 = F_tz1(t(j),z(j));
        k2 = F_tz1(t(j+1),z(j)+dt*k1);

        z(j+1) = z(j) + dt/2*(k1+k2);
    end

    for j = 1:(length(t)-1)
        k1 = F_tz2(t(j),z(j));
        k2 = F_tz2(t(j+1),z(j)+dt*k1);

        z(j+1) = z(j) + dt/2*(k1+k2);
    end

    error(2,k) = abs(z(1)-z(end));

    figure(1)
    loglog(error(2,k),n,'square','color',C(k,:))
    hold on

    %% RK4 (blue)
    % y(t_i+1) = y(t_i) + dt/6*( k1 + 2*k2 + 2*k3 + k4)

    for j = 1:(length(t)-1)
        k1 = F_tz1(t(j),z(j));
        k2 = F_tz1(t(j)+dt/2,z(j)+dt/2*k1);
        k3 = F_tz1(t(j)+dt/2,z(j)+dt/2*k2);
        k4 = F_tz1(t(j)+dt,z(j)+dt*k3);

        z(j+1) = z(j) + dt/6*(k1+2*k2+2*k3+k4);
    end

    for j = 1:(length(t)-1)
        k1 = F_tz2(t(j),z(j));
        k2 = F_tz2(t(j)+dt/2,z(j)+dt/2*k1);
        k3 = F_tz2(t(j)+dt/2,z(j)+dt/2*k2);
        k4 = F_tz2(t(j)+dt,z(j)+dt*k3);

        z(j+1) = z(j) + dt/6*(k1+2*k2+2*k3+k4);
    end

    error(3,k) = abs(z(1)-z(end));

    figure(1)
    loglog(error(3,k),n,'pentagram','color',C(k,:))
    hold on
    
    k = k+1;
end

figure(1)
legend({'FEuler','RK2','RK4'},'Location','best')
xlabel('Error')
ylabel('Time Steps')
title('Numerical Solution Errors for Various Time Steps')
