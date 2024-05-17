% SF1546 - VT24 
% Laboration 1 
% Nikolaos Timoudas % Erik Simert

close all
clear all
clc

%% Uppgift 3 - Samma icke-linjära skalära ekvation med sekantmetoden

% a) 

%Startvärden 
x1a = -1.5; 
x1b = -1;

x2a = -0.3;
x2b = -0.29;

x3a = 0.1;
x3b = 0.2;

x4a = 5;
x4b = 7;

[r1, rel_error1, error1] = sekant(x1a, x1b);
[r2, rel_error2, error2] = sekant(x2a, x2b);
[r3, rel_error3, error3] = sekant(x3a, x3b);
[r4, rel_error4, error4] = sekant(x4a, x4b);

startpunkter1 = [x1a; x2a; x3a; x4a];
startpunkter2 = [x1b; x2b; x3b; x4b];
slutpunkter = [r1(end); r2(end); r3(end); r4(end)];
T1 = table(startpunkter1, startpunkter2, slutpunkter, 'Variablenames', {'Startvärden 1', 'Startvärden 2', 'Rötter'});
disp(T1)

%% b) Skillnad med rötterna från Newtons metod

% Newtons metod hittade inte x1  
% Sekantmetoden hittade inte x2 

%% c) Konvergens - definition

% sida 61 Sauer 
% error(i+1) = |f''(r) / 2f'(r)|*error(i)*error(i-1)

% Tumregel - mellan linjär och kvadratisk konvergens

% K = error(i) / (error(i-1)*error(i-2))

%% d) Konvergens-konstant för x4

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

%% e) Metod som föredras 

% Sekantmetoden föredras i detta tillfälle då 
% derivatan är svår att ta fram

%% Sekantmetoden 

function y = f(x)
y = 61.*x-((x.^2+x+0.03)./(3.*x+1)).^7-20.*x.*exp(-x);
end

function [r_list, rel_error, error] = sekant(x_start0, x_start1)
    
    x0 = x_start0;
    x1 = x_start1;
    tol_error = 1e-8; % tolerans
    n = 100; % iterationer 
    rel_error = zeros(n, 1);
    error = zeros(n, 1);
    r_list = zeros(n, 1);

    for i = 1:n
        f_x0 = f(x0);
        f_x1 = f(x1);
        x_n = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);
        
        r_list(i) = x_n;        
        rel_error(i) = abs((x_n - x1) / x_n) * 100; % relativt fel
        error(i) = abs(x_n - x1);                   % fel
        
        if abs(rel_error(i)) < tol_error
            r_list = r_list(1:i);
            rel_error = rel_error(1:i); 
            error = error(1:i);
            break;
        else
            x0 = x1;
            x1 = x_n;
        end
    end
end
