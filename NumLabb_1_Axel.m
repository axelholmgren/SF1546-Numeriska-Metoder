%% Uppg 2 Icke-linjär skalär ekvation med Newtons metod
% Definierar x värden
x = -3:1e-4:6.5;

% Funktion
f = @(x) 61*x - ((x.^2 + x + 0.03)./(3 * x + 1 )).^7 - 20*x.*exp(-x);
y = f(x);

% Plottar functionen
plot(x, y, 'LineWidth', 2);
ylim([-100, 300])
xlim([-1,11])
title('Function Lab1 Uppgift 2')
xlabel('x');
ylabel('y');
grid on;

% Visa plott
legend('Lab1uppg2func');