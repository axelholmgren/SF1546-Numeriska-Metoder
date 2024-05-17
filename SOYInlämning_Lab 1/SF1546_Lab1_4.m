% SF1546 - VT24 
% Laboration 1 
% Nikolaos Timoudas % Erik Simert

close all
clear all
clc
format long 

%% Uppgift 4 - Interpolation och linjära minsta kvadratmetoden

t1 = [datetime(2024,01,01);
    datetime(2024,02,01);
    datetime(2024,03,01);
    datetime(2024,04,01);
    datetime(2024,05,01);
    datetime(2024,06,01);
    datetime(2024,07,01);
    datetime(2024,08,01);
    datetime(2024,09,01);
    datetime(2024,10,01);
    datetime(2024,11,01);
    datetime(2024,12,01)
    datetime(2024,12,31)];

x = day(t1, 'dayofyear');

time = [6 15; 
        8 06; 
        10 32; 
        13 15; 
        15 55; 
        18 03; 
        18 24; 
        16 38; 
        14 04; 
        11 24; 
        8 46; 
        6 36; 
        6 14];

 y = time(:,1) * 60 + time(:,2);
 
 %% a) 

p_a = polyfit(x, y, length(x)-1); %polynomets grad är antal punkter -1
x_a = linspace(min(x), max(x));
y_a = polyval(p_a, x_a);

subplot(4,2,1);
plot_a = plot(x, y, 'o', x_a, y_a);
xlim([min(x), max(x)])
title('Plot A - Interpolationspolynom');
ylabel('Antal minuter sol');
xlabel('Dag på året');

%% b)

x_b = linspace(min(x),max(x));
y_b = interp1(x,y,x_b,'linear');

subplot(4,2,2);
plot_b = plot(x, y, 'o', x_b, y_b);
xlim([min(x), max(x)])
title('Plot B - Styckvis linjär interpolation');
ylabel('Antal minuter sol');
xlabel('Dag på året');

%% c)

sp_c = spline(x,y);
x_c = linspace(min(x), max(x));
y_c = ppval(sp_c, x_c);

y_c_max = max(y_c)
[y_c_max, max_index_c] = max(y_c);
x_c_max = round(x_c(max_index_c))
date_c = datetime(2024, 1, 1) + days(x_c_max-1)

subplot(4,2,3)
plot_c = plot(x, y, 'o', x_c, y_c);
xlim([min(x), max(x)])
title('Plot C - Splines-approximation');
ylabel('Antal minuter sol');
xlabel('Dag på året');

%% d)

x_d = x(6:8);   % juni till august
y_d = y(6:8);   % juni till august  
p_d = polyfit(x_d, y_d, 2);
y_d_approx = polyval(p_d, x);

subplot(4,2,4)
plot_d = plot(x_d, y_d, 'o', x, y_d_approx);
xlim([min(x), max(x)])
title('Plot D - Andragradspolynom');
ylabel('Antal minuter sol');
xlabel('Dag på året');

%% e)

x_e = x(4:9);    % april till sep
y_e = y(4:9);    % april till sep

p_e = polyfit(x_e, y_e, length(x_e)-1);
y_e_approx = polyval(p_e, x);

subplot(4,2,5)
plot_e = plot(x_e, y_e, 'o', x, y_e_approx);
xlim([min(x), max(x)])
title('Plot E - Minstakvadratmetoden');
ylabel('Antal minuter sol');
xlabel('Dag på året');

%% f)

p_f = polyfit(x, y, length(x)-1);
y_f_approx = polyval(p_f, x);

subplot(4,2,6)
plot_f = plot(x, y_f_approx, 'o', x, y);
xlim([min(x), max(x)])
title('Plot F - Minstakvadratmetoden');
ylabel('Antal minuter sol');
xlabel('Dag på året');

%% g)

f = @(c, x) c(1) + c(2)*cos(2*pi*x/365) + c(3)*sin(2*pi*x/365);
c0 = [1, 1, 1];
c = lsqcurvefit(f, c0, x, y);
y_g_approx = f(c, x);

y_g_max = max(y_g_approx)
[y_g_max, max_index_g] = max(y_g_approx);
x_g_max = x(max_index_g);
date_c = datetime(2024, 1, 1) + days(x_g_max-1)


subplot(4,2,[7,8])
plot(x, y, 'o', x, y_g_approx);
xlim([min(x), max(x)])
title('Plot G - Minstakvadratmetoden');
ylabel('Antal minuter sol');
xlabel('Dag på året');

%% Svar)

%a) För ansats A) krävdes det ett polynom av n:te graden vilket 
% ger n+1 koefficienter, alltså 14. 

%b) De fyra ansatserna är D), E), F), G). 
% D, E, F är ett andragradspolynom y = c1 + c2*x+c3*x^2
% G är givet som y = c1 + c2*cos(wx)+c3*sin(wx)

%c) Metod C) är nog bäst pga använder data från hela året och splines 
%använder kubisk polynom. Soltiden blir max(y_c) = 1111,4 minuter

%d) Metod G) är nog bäst pga periodtiden w = 2pi/365 som tar hänsyn för
%alla årets dagar 

%e) 

%f) M