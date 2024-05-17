% SF1546 - VT24 
% Laboration 2 
% Nikolaos Timoudas % Erik Simert

clear 
close all
clc 
format long

%% 5. Inskjutningsmetoden på samma problem

% Längd 
L = 3.60;   
L0 = 0;

% Randvillkor 
T0 = 310;
TL = 450;

% Värmeledningsförmåga
k = @(x) 3+x./7;
k_prim = 1/7;

% Värmemängd 
Q = @(x) 280.*exp(-(x-L./2).^2);

% ODE-system som ska lösas 
ode_system = @(x, u) [u(2); -1./k(x) .* (Q(x) + k_prim .* u(2)) ];

% Inskjutningsmetoden
ode_options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
shooting_function = @(T_guess) ode45(@(x, u) ode_system(x, u), [0, L], [T0, T_guess], ode_options).y(1, end) - TL;

% Gissning 
T_guess = fzero(shooting_function, TL);

% Lösning med gissningen 
[x, T] = ode45(@(x, u) ode_system(x, u), [0, L], [T0, T_guess], ode_options);

% T vid x = 1.65
[~, x_temp] = min(abs(x - 1.65));
T_wanted = T(x_temp);

T1 = table(T_wanted, 'VariableNames', {'T vid x = 1.65'});
disp(T1)

% Plot
figure;
plot(x, T(:, 1), 'LineWidth', 1.3);
xlabel('Position [m]');
ylabel('Temperatur [K]');
title('Temperaturfördelning');
xlim([0 L])
grid on;
set(gca,'FontSize',16);
set(gca,'FontName','times');
