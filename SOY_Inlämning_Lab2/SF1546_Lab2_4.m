% SF1546 - VT24 
% Laboration 2 
% Nikolaos Timoudas % Erik Simert

clear 
close all
clc 
format long

%% 4. Differentialekvationer - Randvärdesproblem 

% b) 

% Längd 
L = 3.60;   
L0 = 0;

% Randvillkor 
T0 = 310;
TL = 450;

% Värmeledningsförmåga
k = @(x) 3+x./7;
% Värmemängd 
Q = @(x) 280.*exp(-(x-L./2).^2);

% Inre punkter 
n = [1e2, 1e3, 1e4, 1e5, 1e6];

T_cell = cell(size(n));
x_cell = cell(size(n));
T_max = zeros(size(n));
T_min = zeros(size(n));
T_medel = zeros(size(n));
T_wanted = zeros(size(n));
error = zeros(size(n));
for i = 1:length(n)
    [T, x, D] = T_calc(L0, L, T0, TL, n(i), k, Q);
    T_cell{i} = T;
    x_cell{i} = x;
    T_max(i) = max(T);
    T_min(i) = min(T);
    T_medel(i) = mean(T);

    [~, idx] = min(abs(x - 1.65));
    T_wanted(i) = T(idx);

    % Beräkna felgränsen 
    if i > 1
        error(i) = abs(T_wanted(i) - T_wanted(i-1));
    end
end

% Tabell
T1 = table(n', T_wanted', error', ...
     'VariableNames', {'Inre punkter', 'T när x = 1.65', 'Felgräns'});
disp(T1);

% Plot 
figure(1); 
plot(x, [T(1);T;T(end)], 'LineWidth', 1.3);
grid on;
xlabel('Avstånd [m]');
ylabel('Temperatur [K]');
xlim([0 L]);
title('Temperaturfördelning');
set(gca,'FontSize',16);
set(gca,'FontName','times');

%% Funktion som beräknar temperaturfördelningen 

function [T, x, D] = T_calc(L0, L, T0, TL, n, k, Q)

    % Diskretisering
    h = (L-L0)./(n+1);

    % Steg vektor och inre punkter 
    x = linspace(L0,L,n+2);
    x_inre = x(2:end-1)'; 
    k_values = k(x_inre);
    
    % A matrisen 
    diag_upper = ((-1 ./ 14*h) - (k_values ./ h.^2)); 
    diag_center = ((2 .* k_values) ./ h.^2);
    diag_under = ((1 ./ 14*h) - (k_values ./ h.^2)); 
    
    D = spdiags([diag_under, diag_center, diag_upper], -1:1, n, n);
    D(1, :) = 0;
    D(1, 1) = 1;
    D(end, :) = 0;
    D(end, end) = 1;

    b = Q(x_inre);
    b(1) = T0;
    b(end) = TL;
    
    % Löser ut T 
    T = D \ b;
end