% Zane Grothe
% AERO 7970
% ODE Solver
% 10/14/22

clear all
close all
clc

%   Goal
%        dy/dt = f(t,y)
%
%   dy/dt = cos(t)
%   y(t) = sin(t) + K
%   y(t0) = sin(t0) + K
%   y(t) = sin(t) + y(t0) - sin(t0)

%  Initial conditions
t0 = 0;                  % starting time
tf = 25;                 % final time
dt = 0.05;               % step size
t = t0:dt:tf;            % range of result
y = zeros(1,length(t));  % preallocation
y(1) = 1;                % initial condition
F_ty = @(t,y) cos(t);    % dy/dt

% Color matrix to pick from for plotting
C=[1,0,0; 0,1,0; 0,0,1; .929,.694,.125; 0,1,1; 1,0,1; 0,0,0; .85,.325,.098];
% [ red ; green; blue ; gold          ; cyan ; mag. ; black; brown        ]

%% Exact (red)
% y(t) = y(t0) + sin(t) - sin(t0)

y(1) = F_ty(t(1),y(1));
for i = 1:(length(t)-1)
    y(i+1) = y(1) + sin(t(i+1)) - sin(t(1));
end

y_exact = y;

figure(1)
plot(t,y_exact,'color',C(1,:))
hold on

%% Forward Euler (green)
% y(t_i+1) = y(t_i) + dt*f(t_i, y(t_i))

y(1) = F_ty(t(1),y(1));
for i = 1:(length(t)-1)
    y(i+1) = y(i) + dt*F_ty(t(i),y(i));
end

y_feuler = y;

figure(1)
plot(t,y_feuler,'color',C(2,:))
hold on

%% Backward Euler (blue)
% y(t_i+1) = y(t_i) + dt*f(t_i+1, y(t_i+1))

y(1) = F_ty(t(1),y(1));
for i = 1:(length(t)-1)
    tol = 1e-9;         % tolerance
    max_iter = 1000;    % maximum number of picard iterations
    yg = 1;             % y is a GUESS
    iter = 0;
    y_old = yg;
    error = 1.0;        % Let's set the initial error huge
    
    while error > tol && iter < max_iter
        y_new = yg + dt*F_ty(t(i+1), y_old);
        
        error = abs(y_new - y_old) / abs(y_old);
        y_old = y_new;
        
        iter = iter + 1;
    end
    
    if error > tol
        disp('Picard iteration has failed')
    end
    
    y(i+1) = y(i) + dt*F_ty(t(i+1),y_new);
end

y_beuler = y;

figure(1)
plot(t,y_beuler,'color',C(3,:))
hold on

%% Crank Nicholson (gold)
% y(t_i+1) = y(t_i) + dt/2*( f(t_i,y(t_i)) + f(t_i+1, y(t_i+1)) )

y(1) = F_ty(t(1),y(1));
for i = 1:(length(t)-1)
    tol      = 1e-9;    % tolerance
    max_iter = 1000;    % maximum number of picard iterations
    yg       = 1;       % initial y (guess)
    iter     = 0;
    y_old    = yg;
    error    = 1.0;     % initial error (large)
    
    while error > tol && iter < max_iter
        y_new = yg + dt/2*(F_ty(t(i), y(i))+F_ty(t(i),y_old));
    
        error = abs(y_new - y_old) / abs(y_old);
        y_old = y_new;
    
        iter = iter + 1;
    end

    if error > tol
        disp('Picard iteration has failed')
    end
    
    y(i+1) = y(i) + dt/2*(F_ty(t(i),y(i))+F_ty(t(i+1),y_new));
end

y_crank = y;

figure(1)
plot(t,y_crank,'color',C(4,:))
hold on

%% RK2 (cyan)
% y(t_i+1) = y(t_i) + dt/2*( k1 + k2 )

y(1) = F_ty(t(1),y(1));
for i = 1:(length(t)-1)
    k1 = F_ty(t(i),y(i));
    k2 = F_ty(t(i+1),y(i)+dt*k1);

    y(i+1) = y(i) + dt/2*(k1+k2);
end

y_rk2 = y;

figure(1)
plot(t,y_rk2,'color',C(5,:))
hold on

%% RK4 (magenta)
% y(t_i+1) = y(t_i) + dt/6*( k1 + 2*k2 + 2*k3 + k4)

y(1) = F_ty(t(1),y(1));
for i = 1:(length(t)-1)
    k1 = F_ty(t(i),y(i));
    k2 = F_ty(t(i)+dt/2,y(i)+dt/2*k1);
    k3 = F_ty(t(i)+dt/2,y(i)+dt/2*k2);
    k4 = F_ty(t(i)+dt,y(i)+dt*k3);
    
    y(i+1) = y(i) + dt/6*(k1+2*k2+2*k3+k4);
end

y_rk4 = y;

figure(1)
plot(t,y_rk4,'color',C(6,:))
hold on
xlabel('t')
ylabel('y')
title('Numerical Solutions for dy/dt=cos(t)')
legend({'Exact','FEuler','BEuler','Crank','RK2','RK4'},'Location','northeast')
