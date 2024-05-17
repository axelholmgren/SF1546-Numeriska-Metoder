% SF1546 - VT24 
% Laboration 1 
% Nikolaos Timoudas % Erik Simert

close all
clear 
clc

format long 

%% Uppgift 6 - Integral med förbehandling

%a) 

x_a = linspace(0, 1e-4, 100);
f_a = @(x_a) (1 - exp(-(x_a/5).^3))./(5*x_a.^3);

subplot(2,1,1)
plot(x_a, f_a(x_a))
xlabel('x')
ylabel('Integranden')
title('A - Integranden nära 0')

%% b) 

x_b = linspace(1e-4, 1000, 1000);
f_b =@(x_b) (1 - exp(-(x_b/5).^3))./(5*x_b.^3);

subplot(2,1,2)
plot(x_b, f_b(x_b))
xlabel('x')
ylabel('Integranden')
title('B - Svanskapning')

%% c)

%Gränsvärde mot noll - konstant som läggs som undre gräns 

f_zero = 1/625;

%% d)

%Det krävs inga ytterligare förberedelser för att beräkna integralen. 
%Trapetsregeln kräver att integralen är jämn 

%% e) 

%Trapetsregeln
lower_bound = f_zero;
upper_bound = 1000;
step_size = 1e6;
x = linspace(lower_bound, upper_bound, step_size);
y = f_b(x);
dx = x(2) - x(1);

%Arean av varje trapets sen summera totala arean
trapets_area = (y(1:end-1) + y(2:end)) / 2 * dx;
f_int_trapets = sum(trapets_area);

%% f)

h = 1e6;
steg_n = [h, h/2, h/4, h/8, h/16];
f_integral = [];

for i = 1:length(steg_n)
    
    steg = steg_n(i);
    x = linspace(lower_bound, upper_bound, steg);
    y = f_b(x);
    dx = x(2) - x(1);
    trapets_area = (y(1:end-1) + y(2:end)) / 2 * dx;
    f_int_trapets = sum(trapets_area);
    
    f_integral(i) = f_int_trapets;
end

fprintf('\nSvar:\n')
fprintf('%-12s %-12s\n', 'Steglängd', 'Integral');
for i = 1:length(f_integral)
    fprintf('%-12.0f %-12.16f\n', steg_n(i), f_integral(i));
end

fprintf('\nFelgränser [x10^-12]: \n');
err1 = (f_integral(1)-f_integral(2));
err2 = (f_integral(2)-f_integral(3));
err3 = (f_integral(3)-f_integral(4));
err4 = (f_integral(4)-f_integral(5));
error_tab = abs([err1, err2, err3, err4]*10^12);
header = {'h-h/2', 'h/2-h/4', 'h/4-h/8', 'h/8-h/16'};
fprintf('%12s %12s %12s %12s\n', header{:});
fprintf('%12.6f %12.6f %12.6f %12.6f\n', error_tab);