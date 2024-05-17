% SF1546 - VT24 
% Laboration 2 
% Nikolaos Timoudas % Erik Simert

clear 
close all
clc 

format long 
%% 3. Differentialekvationer - Begynnelsevärdesproblem 

%a) 

% Andra ordningens differentialekvation 
% φ'' + (g/L)*sin(φ) = 0 

%System av första ordningens differentialekvationer 
% φ' = w
% w' = -(g/L)*sin(φ)

%Begynnelsevärden 
phi_0 = (6*pi)/7;
phi_prim_0 = 0.8;

%% b) 

%Konstanter
T = 20;
L = 2.5;
g = 9.81;

% φ motsvarar y(1) (vinkel)
% φ' motsvarar y(2) (vinkelhastighet)

%System av ODE 
ode_system = @(t, y) [y(2); -(g/L)*sin(y(1))];

%Initial Conditions
IC = [phi_0; phi_prim_0];  

%ode45 lösning
[t, y] = ode45(ode_system, [0, T], IC);
phi = y(:, 1); %Vinkel - första kolumnen
phi_prim = y(:, 2); %Vinkelhastighet - andra kolumnen

% Plotta vinkel och vinkelhastighet
figure(1)
subplot(2, 1, 1);
p1 = plot(t, phi, 'LineWidth', 1.3);
xlabel('Tid (s)');
ylabel('Vinkel (rad)');
title('Vinkel som funktion av tid');
set(gca,'FontSize',16);
set(gca,'FontName','times');

subplot(2, 1, 2);
p2 = plot(t, phi_prim, 'LineWidth', 1.3);
xlabel('Tid (s)');
ylabel('Vinkelhastighet (rad/s)');
title('Vinkelhastighet som funktion av tid');
set(gca,'FontSize',16);
set(gca,'FontName','times');

%% d) 

[T1] = period_calc(t, phi);

%% e) 
%Desto längre pendel --> längre period 
% T = 2pi*sqrt(L/g) 

%Konstanter
L_new = 2.7;

%System av ODE 
ode_system_new = @(t, y_new) [y_new(2); -(g/L_new)*sin(y_new(1))];

%ode45 lösning
[t, y_new] = ode45(ode_system_new, [0, T], IC);
phi_new = y_new(:, 1); %Vinkel - första kolumnen
phi_prim_new = y_new(:, 2); %Vinkelhastighet - andra kolumnen

[T2] = period_calc(t, phi_new);

T1 = table(T1, T2, 'VariableNames', {'Period med L = 2.5', 'Period med L = 2.7'});
disp(T1)

% c) 
figure(2)
anim(t, phi, L); 

%% Beräkning av periodtiden

function [T_period] = period_calc(t, phi)
    %Interpolera svängningen 
    t_spline = linspace(t(1), t(end), 1e6); % t vektor 
    phi_spline = spline(t, phi, t_spline); % Interpolation
    
    [~, indices] = sort(abs(phi_spline));
    phi_zero = phi_spline(indices(1:4));
    t_find = t_spline(indices(1:4));
    t_sort = sort(t_find);

    % Beräkna svingningstiden 
    T_diff = diff(t_sort); % Skillnaden mellan tider ger perioden
    T_period = 2 * mean(T_diff); % Beräkna den genomsnittliga perioden och multiplicera med 2 för att få svingningstiden
end

%% Animering av pendeln  

function anim(tut,fiut,L)
    figure(2)
    for i=1:length(tut)-1
        x0=L*sin(fiut(i));y0=-L*cos(fiut(i));
        plot([0,x0],[0,y0],'-o', 'LineWidth', 1.3)
        title('Animering av pendeln')
        set(gca,'FontSize',16);
        set(gca,'FontName','times');
        axis('equal')
        axis([-1 1 -1 0]*1.2*L)
        drawnow
        pause(tut(i+1)-tut(i))
    end
end

