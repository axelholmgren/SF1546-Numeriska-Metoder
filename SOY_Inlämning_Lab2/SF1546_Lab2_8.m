% SF1546 - VT24 
% Laboration 2 
% Nikolaos Timoudas % Erik Simert

close all 
clear 
clc

%% 8. Mera kinematik

% a) 

% Parametrar
s0 = 2.2; 
e0 = 2.8;
c0 = 0;
p0 = 0;
k1 = 11.0;
k2 = 1.6;

% Time span
t_step = 25;
t_vec = linspace(0, 1.5, t_step);

% Initial conditions
IC = [s0; e0; c0; p0];

steg = [1e-2, 1e-3, 1e-4, 1e-5];
k3_start = 1;
k3_end = 1.5;
p_target = 1.75;

for i=1:length(steg)
    n = round( (k3_end - k3_start) / steg(i));
    k3_vec = linspace(k3_start, k3_end, n);
    p_end = size(k3_vec);


    for j=1:length(k3_vec)
        [t, y] = ode45(@(t, y) ode_system(t, y, k1, k2, k3_vec(j)), t_vec, IC);
        p_end(j) = y(end, 4);

        if p_end(j) >= p_target
            k3_index = j;
            k3_new = k3_vec(k3_index);
            break;
        end
    end
    
    resultat(i, 1) = steg(i);
    resultat(i, 2) = p_end(end);
    resultat(i, 3) = k3_new;
end

error1_k3= resultat(1,3) - resultat(2,3);
error2_k3= resultat(2,3) - resultat(3,3);
error3_k3= resultat(3,3) - resultat(4,3);
error_k3 = [0 ; error1_k3; error2_k3; error3_k3];

T2 = table(resultat(:,1), resultat(:,2), resultat(:,3), error_k3, 'VariableNames', {'Steg', 'Koncentration', 'k3', 'Felgräns'});
disp(T2)

koncentration = {'s'; 'e'; 'c'; 'p'};
y_end = [y(end, 1); y(end, 2); y(end, 3); y(end, 4)];
T1 = table(koncentration, y_end, 'VariableNames', {'Ämne', 'Koncentration'});
disp(T1)

% Plot 
figure(1); hold on; grid on;
plot(t, y(:, 1), 'LineWidth', 1.3);
plot(t, y(:, 2), 'LineWidth', 1.3);
plot(t, y(:, 3), 'LineWidth', 1.3);
plot(t, y(:, 4), 'LineWidth', 1.3);
xlabel('Tid'); ylabel('Koncentration');
title('System');
legend('s', 'e', 'c', 'p');
hold off

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
