% SF1546 - VT24 
% Laboration 2 
% Nikolaos Timoudas % Erik Simert

close all 
clear all
clc

%% 7. Kemisk dynamik - kinetik

% b) 

% Parametrar
s0_b = 2.9; 
e0_b = 1.3;
c0_b = 0;
p0_b = 0;
k1_b = 1.0;
k2_b = 0.8;
k3_b = 1.1;

% Tidsvektor 
dt = 100;
t_vec = linspace(0, 1.5, dt);

% Initial conditions
IC_b = [s0_b; e0_b; c0_b; p0_b];

% Löser systemet 
[t, y_b] = ode45(@(t, y_b) ode_system(t, y_b, k1_b, k2_b, k3_b), t_vec, IC_b);

% Plot 
subplot(2,2,1); hold on; grid on;
plot(t, y_b(:, 1), 'LineWidth', 1.3);
plot(t, y_b(:, 2), 'LineWidth', 1.3);
plot(t, y_b(:, 3), 'LineWidth', 1.3);
plot(t, y_b(:, 4), 'LineWidth', 1.3);
xlabel('Tid'); ylabel('Koncentration');
title('System B');
legend('s', 'e', 'c', 'p');
hold off

%% c) 

% Parametrar
s0_c = 0.4; 
e0_c = 0.6;
c0_c = 0;
p0_c = 0;
k1_c = 2.1;
k2_c = 1.2;
k3_c = 1.5;

% Initial conditions
IC_c = [s0_c; e0_c; c0_c; p0_c];

% Löser systemet 
[t, y_c] = ode45(@(t, y_c) ode_system(t, y_c, k1_c, k2_c, k3_c), t_vec, IC_c);

% Plot 
subplot(2,2,2); hold on; grid on;
plot(t, y_c(:, 1), 'LineWidth', 1.3);
plot(t, y_c(:, 2), 'LineWidth', 1.3);
plot(t, y_c(:, 3), 'LineWidth', 1.3);
plot(t, y_c(:, 4), 'LineWidth', 1.3);
xlabel('Tid'); ylabel('Koncentration');
title('System C');
legend('s', 'e', 'c', 'p');
hold off

%% d)

% Parametrar
s0_d = 2.1; 
e0_d = 2.9;
c0_d = 0;
p0_d = 0;
k1_d = 10.0;
k2_d = 1.6;
k3_d = 1.3;

% Initial conditions
IC_d = [s0_d; e0_d; c0_d; p0_d];

% Löser systemet 
[t, y_d] = ode45(@(t, y_d) ode_system(t, y_d, k1_d, k2_d, k3_d), t_vec, IC_d);

% Plot 
subplot(2,2,3); hold on; grid on;
plot(t, y_d(:, 1), 'LineWidth', 1.3);
plot(t, y_d(:, 2), 'LineWidth', 1.3);
plot(t, y_d(:, 3), 'LineWidth', 1.3);
plot(t, y_d(:, 4), 'LineWidth', 1.3);
xlabel('Tid'); ylabel('Koncentration');
title('System D');
legend('s', 'e', 'c', 'p');
hold off

%% e)

% Parametrar
s0_e = 1.1; 
e0_e = 2.1;
c0_e = 0;
p0_e = 0;
k1_e = 1.8;
k2_e= 1.0;
k3_e = 0.02;

% Initial conditions
IC_e = [s0_e; e0_e; c0_e; p0_e];

% Löser systemet 
[t, y_e] = ode45(@(t, y_e) ode_system(t, y_e, k1_e, k2_e, k3_e), t_vec, IC_e);

% Plot 
subplot(2,2,4); hold on; grid on;
plot(t, y_e(:, 1), 'LineWidth', 1.3);
plot(t, y_e(:, 2), 'LineWidth', 1.3);
plot(t, y_e(:, 3), 'LineWidth', 1.3);
plot(t, y_e(:, 4), 'LineWidth', 1.3);
xlabel('Tid'); ylabel('Koncentration');
title('System E');
legend('s', 'e', 'c', 'p');
hold off

koncentration = {'s'; 'e'; 'c'; 'p'};
y_b_end = [y_b(end, 1); y_b(end, 2); y_b(end, 3); y_b(end, 4)];
y_c_end = [y_c(end, 1); y_c(end, 2); y_c(end, 3); y_c(end, 4)];
y_d_end = [y_d(end, 1); y_d(end, 2); y_d(end, 3); y_d(end, 4)];
y_e_end = [y_e(end, 1); y_e(end, 2); y_e(end, 3); y_e(end, 4)];
T1 = table(koncentration, y_b_end, y_c_end, y_d_end, y_e_end, 'VariableNames', {'Ämne', 'b)', 'c)', 'd)', 'e)'});
disp(T1)

%% ODE system 

function dydt = ode_system(t, y, k1, k2, k3)
    s = y(1);
    e = y(2);
    c = y(3);
    
    dsdt = -k1 * s * e + k2 * c;
    dedt = -k1 * s * e + k2 * c + k3 * c;
    dcdt = k1 * s * e - k2 * c - k3 * c;
    dpdt = k3 * c;
    
    dydt = [dsdt; dedt; dcdt; dpdt];
end
