% SF1546 - VT24 
% Laboration 2 
% Nikolaos Timoudas % Erik Simert

clear 
close all
clc 

%% 2. Numerisk integration: Rotationssymmetrisk lur

% a) & b) 

% Steglängder
h = 0.5;
steg_h = [h ; h/2 ; h/4 ; h/8 ; h/16 ; h/32 ; h/64; h/128 ; h/256 ; h/512 ; h/1024 ; h/2056];

% Längd
L = 4; 

% Beräkna volym 
[V, error, ~, ~] = volume_calc(L, steg_h);

%% c) 

% Ny volym
p = 0.72;
V_target = V(end) * p;

% Startvärde för L 
L_0 = 3;
L_1 = 1;
% Tolerans 
tol = 1e-6;

L_new_temp = size(steg_h);
L_error = size(steg_h);
L_error_temp = size(steg_h);

for i = 3:length(steg_h)

    [L_new] = calc_new_L(tol, V_target, L_0, L_1, steg_h(i));
    L_new_temp(i) = L_new;
    L_new_values = L_new_temp';
    L_error_temp(i) = abs(L_new_values(i)-L_new_values(i-1));
    L_error = L_error_temp';
end

% Tabell över resultat  
T1 = table(steg_h, V', error', 'VariableNames', {'Steglängd', 'Volym', 'Felgräns'});
disp(T1)

T2 = table(steg_h(4:end,:), L_new_values(4:end,:), L_error(4:end,:), 'VariableNames', {'Steglängd', 'Längd för 0.72*V', 'Felgräns'});
disp(T2)

%% d) 

[~, ~, x_vec, f_vec] = volume_calc(L, steg_h);

n = 50;
x_points = round(linspace(1, length(x_vec), n));
x = x_vec(x_points);
f = f_vec(x_points);
        
phi = linspace(0, 2*pi, length(x));
X = x' * ones(size(phi)); 
Y = f' * cos(phi); 
Z = (f' * sin(phi))';  

figure(1);
surf(X, Y, Z);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D-bild av luren');

%% Funktion som hittar nya L

function [L_new] = calc_new_L(tolerans, V_target, L_0, L_1, steg_h)

    V_fun = @(L) volume_calc(L, steg_h);
    
    while abs(L_1 - L_0) > tolerans
        
      % Beräkna funktionsvärdena för startvärdena
      f_0 = V_fun(L_0) - V_target;
      f_1 = V_fun(L_1) - V_target;

      % Beräkna det nya värdet för L
      L_new = L_1 - f_1 * (L_1 - L_0) / (f_1 - f_0);
      % Uppdatera startvärdena
      L_0 = L_1;
      L_1 = L_new;
    end
end    

%% Funktion som beräknar volym 
    
function [V, error, xVec, y] = volume_calc(L, steg_h)

    % Differentialekvation
    f = @(x, y) -(1/6 + (pi * sin(pi * x)) / (1.6 - cos(pi * x))) .* y;
    x0 = 0;
    x1 = L;

    for i = 1:length(steg_h)

        n = round((x1 - x0) / steg_h(i));

        xVec = linspace(x0, x1, n+1);
        
        % Randvärde 
        y(1) = 2.5;

        for j = 1:n
            y(1 + j) = y(j) + steg_h(i) * f(xVec(j), y(j)); % Euler method
        end

        % Trapetsregeln
        g = (y.^2);
        V(i) = pi .* trapz(xVec, g);

        % Felgräns 
        if i > 1
            error(i) = abs(V(i) - V(i-1));
        end
    end
end
