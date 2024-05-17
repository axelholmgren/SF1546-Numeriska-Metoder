% SF1546 - VT24 
% Laboration 1 
% Nikolaos Timoudas % Erik Simert

close all
clear 
clc

format long 

%% Uppgift 7 - Säker integral

%a)

x = linspace(0,6, 1000);
f = @(x) 153*exp(-((11*x-pi)./0.004).^2);
y = f(x);
tol = 1e-12;
plot(x, y) 
xlabel('x')
ylabel('y')
set(gca,'FontSize',16);
set(gca,'FontName','times');

%% Quad 

[ans_quad1] = calc_quad(f, 0, 0.28, tol);
[ans_quad2] = calc_quad(f, 0.28, 0.3, tol);
[ans_quad3] = calc_quad(f, 0.3, 6, tol);

ans_quad = (ans_quad1 + ans_quad2 + ans_quad3);

%% Integral

[ans_integral1] = calc_integral(f, 0, 0.285, tol);
[ans_integral2] = calc_integral(f, 0.285, 0.29, tol);
[ans_integral3] = calc_integral(f, 0.3, 6, tol);

ans_integral = (ans_integral1 + ans_integral2 + ans_integral3);

%% Felgräns

E_quad = tol*abs(ans_quad);
E_integral = tol*abs(ans_integral);

answers = [ans_quad ; ans_integral];
errors = [E_quad ; E_integral];
type = ["quad"; "integral"];

T1 = table(type, answers, errors, 'VariableNames', {'Metod', 'Svar', 'Felgräns'});
disp(T1)

%% Functions 

function [answer] = calc_quad(f, lower_bound, upper_bound, tol)
    answer = quad(f,lower_bound,upper_bound,tol);
end

function [answer] = calc_integral(f, lower_bound, upper_bound, tol)
    answer = integral(f,lower_bound,upper_bound,'AbsTol', tol);
end

%% b
% Genom att börja med en låg feltolerans och sedan öka denna för att få 
% ett mer exakt svar kan man spara på arbetet. Vidare genom att dela upp
% funktionen där det sker snabba förändringar, exempelvis genom att notera
% ifall funktionen har asymptoter. Dessa områden kan sedan beräknas för sig
% genom att använda någon av följande metoder: substitution, partiell
% integration, uppdelning, kapning. 
% I detta fall delas integralen upp där det finns en asymptot.

