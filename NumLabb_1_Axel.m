%% Uppg 2 Icke-linjär skalär ekvation med Newtons metod
% b 
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

% c
% Newtons Metod
initialgissning = 1;
tolerans = 1e-10;
maxIterationer = 100;

% Definierar derivatan
fDerivata = @(x) 61 + 21*((x.^2+x+0.03)./(3.*x+1)^8) - 7*((2.*x+1)*(x.^2+x+0.03).^6)./((3*x+1).^7) - 20.*exp(-x) + 20*x.*exp(-x);

% Newtons metod
forraApprox = initialgissning;
for iteration = 1:maxIterationer
    x1 = initialgissning - f(initialgissning) / fDerivata(initialgissning);
    if abs(x1 - initialgissning) < tolerans
        break;
    end

    initialgissning = x1;
end

rot = initialgissning;
iterations = iteration;
disp(['Rot: ' num2str(rot)]);
disp(['Antal interationer: ' num2str(iterations)]);
