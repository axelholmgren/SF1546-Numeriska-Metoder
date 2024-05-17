% SF1546 - VT24 
% Laboration 1 
% Nikolaos Timoudas % Erik Simert

close all
clear all
clc
format long 

%% Uppgift 2 - Icke-linjär skalär ekvation med Newtons metod

% a) 

%Använda medelvärdessatsen för att grovlokalisera rötter   

%% b) Plot av funktionen och dess nollställen

x = -2:0.1:7;
fun = 61.*x-((x.^2+x+0.03)./(3.*x+1)).^7-20.*x.*exp(-x);

figure(1)
plot(x,fun);
title('Plot Uppgift 2');
ylabel('f(x)');
xlabel('x');
grid on;


%Vi har rötter enligt wolfram-alpha på:
% xa = -1.1151              x1 = hittar ej 
% xb = -0.2910              x2 = -1 
% xc = 5.3341 * 10^-13      x3 = 0
% xd = 6.3971               x4 = 6

%% c) Newtons metod 

%Startvärden 
x1 = -1;
x2 = -0.3;
x3 = 1;
x4 = 6;

[r1, rel_error1, error1] = newton(x1);
[r2, rel_error2, error2] = newton(x2);
[r3, rel_error3, error3] = newton(x3);
[r4, rel_error4, error4] = newton(x4);

startpunkter = [x1; x2; x3; x4];
slutpunkter = [r1(end); r2(end); r3(end); r4(end)];
T1 = table(startpunkter, slutpunkter, 'Variablenames', {'Startvärden', 'Rötter'});
disp(T1)


%% d) Kvadratisk konvergens - definition 

%  1.24 i Sauer ger:
% e(i+1) = e(i)^2 * |f''(c(i)) / 2 f'(x(i))|
% där c(i) motsvarar ett värde mellan x(i) och roten r 

% Tumregeln
% e(i+1) = M * e(i)^2

%% e) Konvergens-konstanten för x4

E1 = abs(r1(1:end-1)-r1(2:end));
log_E1 = log(E1(1:end-1))-log(E1(2:end));
K1 = (log_E1(2:end) ./ log_E1(1:end-1));

E2 = abs(r2(1:end-1)-r2(2:end));
log_E2 = log(E2(1:end-1))-log(E2(2:end));
K2 = (log_E2(2:end) ./ log_E2(1:end-1));

E3 = abs(r3(1:end-1)-r3(2:end));
log_E3 = log(E3(1:end-1))-log(E3(2:end));
K3 = (log_E3(2:end) ./ log_E3(1:end-1));

E4 = abs(r4(1:end-1)-r4(2:end));
log_E4 = log(E4(1:end-1))-log(E4(2:end));
K4 = (log_E4(2:end) ./ log_E4(1:end-1));

K = [K1(end) ; K2(end) ; K3(end) ; K4(end)];

T2 = table(slutpunkter, K, 'Variablenames', {'Rötter', 'Konvergensordning'});
disp(T2)

%% Implementering av Newtons Metod 

%function y = f(x)
%y = 61.*x-((x.^2+x+0.03)./(3.*x+1)).^7-20.*x.*exp(-x);
%end

%function y = df(x)
%y = ((21.*(x.^2+x+0.03).^7)./(3.*x+1).^8) - ((7.*(2.*x+1).*(x.^2+x+0.03).^6)./(3.*x+1).^7)-20.*exp(-x)+20.*exp(-x)+61;
%end 

function y = f(x)
    e = exp(1);
    y = 61*x - ((x.^2 + x + 0.03)./(3*x + 1 )).^7 - 20*x.*e.^(-x);
end

function y = df(x)
    e = exp(1);
    y = (21.*(x.^2 + x + 0.03).^7)./((3.*x + 1).^8) - (7.*(2.*x + 1).*(x.^2 + x + 0.03).^6)./((3.*x + 1).^7) - 20.*e.^(-x) + (20.*e.^(-x)).*x + 61;

end 

%function [f, d] = f_icke_lin(x)
%    % beräknar värde och derivata av den ickelinjära funktionen i uppg 2 lab 1. returneras som [funktionsvärde , derivata, noll för nollställe]
%    e = exp(1);
%    f = 61*x - ((x.^2 + x + 0.03)./(3*x + 1 )).^7 - 20*x.*e.^(-x);
%    % derivata från wolfram alpha
%    d = (21.*(x.^2 + x + 0.03).^7)./((3.*x + 1).^8) - (7.*(2.*x + 1).*(x.^2 + x + 0.03).^6)./((3.*x + 1).^7) - 20.*e.^(-x) + (20.*e.^(-x)).*x + 61;
%end

function [r_list, rel_error, error] = newton(x_start)
    
    x0 = x_start;
    tol_error = 1e-8;   %tolerans
    n = 100;            %iterationer 
       
    rel_error = zeros(n, 1);
    error = zeros(n, 1);
    r_list = zeros(n, 1);

    for i=1:n
        x_n = x_start - f(x_start)/df(x_start);
        
        r_list(i) = x_n;
        rel_error(i) = abs((x_n-x_start)/(x_n))*100; %relativt fel
        error(i) = abs(x_n - x0); 

        if abs(rel_error(i)) < tol_error
            r_list = r_list(1:i);
            rel_error = rel_error(1:i);
            error = error(1:i);
            break;
        else 
            x_start = x_n;
        end

    end
end


