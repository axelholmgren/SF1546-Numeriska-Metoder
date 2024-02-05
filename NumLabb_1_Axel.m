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

%% uppg 4 Interpolation och linjära minsta kvadratmetoden

% Matrisen A ger datum: tid (sekunder)

A = [1 32 60 91 121 152 182 213 244 274 305 335 365; 374 486 633 794 955 1084 1104 998 845 684 527 396 374]';

% Definierar x och y utifrån matrisen 

x = A(:,1);
y = A(:,2);

% Definierar antal dagar

d = 1:1:365;

% Skapar anpassningar
figure()

% A)
subplot(3,3,1);
    % grad
    grad = length(x) -1;
    
    % Anpassa polynomet till datan
    k1 = polyfit(x,y,grad);
    
    % Skapa polynom med koefficienter
    p1 = polyval(k1,d);
    
    % Plotta
    hold on;
    plot(d,p1);
    plot(x,y,'b.');
    title("Interpolationspolynom genom alla punkter");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0,365,0,1300]);
    hold on;

    
% B)
subplot(3,3,2);
    k2 = interp1(x,y,d,"linear");
    plot(d,k2);
    hold on;
    plot(x,y,'b.');
    title("Styckvis interpolation genom alla punkter");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0,365,0,1300]);
    hold on;


% C)
subplot(3,3,3);
    k3 = spline(x,y,d);
    plot(d,k3);
    hold on;
    plot(x,y,'b.');
    title("Splines-approximation genom alla punkter");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0,365,0,1300]);
    hold on;

    
% D)
subplot(3,3,4);
    delX = x(6:8);
    delY = y(6:8);
    % Grad 2
    k4 = polyfit(delX,delY,2);
    p2 = polyval(k4,d);
    plot(d,p2);
    hold on;
    plot(x,y,'b.');
    title("2a gradspolynom (1 jun-1 aug");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0,365,0,1300]);
    hold on;


% E)
subplot(3,3,5);
    delX2 = x(4:9);
    delY2 = y(4:9);
    % Grad 2
    k5 = polyfit(delX2, delY2, 2);
    p3 = polyval(k5,d);
    plot(d,p3);
    hold on;
    plot(x,y,'b.');
    title("2a gradspolynom (1 apr-1 sept");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0,365,0,1300]);
    hold on;


% F) 
subplot(3,3,6);
    % Grad 2
    k6 = polyfir(x,y,2);
    p4 = polyval(k6,d);
    plot(d,p4);
    hold on;
    plot(x,y,'b.');
    title("2a gradspolynom genom alla punkter");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0,365,0,1300]);
    hold on;


% G)
subplot(3,3,7);
    w= 2*pi/365;
    c1 = ones(13,1);
    c2 = cos(w*x);
    c3 = sin(w*x);

    C = [c1,c2,c3];

