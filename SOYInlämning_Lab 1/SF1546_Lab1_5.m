% SF1546 - VT24 
% Laboration 1 
% Nikolaos Timoudas % Erik Simert

close all
clear 
clc

format long 

%% Uppgift 6 - Tillförlitlighet

x0 = 0.01;
x1 = 6;

[root0] = inc_dec_fun(x0, 1, 1, 1);
[root1] = inc_dec_fun(x1, 1, 1, 1);

%% a) 

[root0_20_inc] = inc_dec_fun(x0, 1.01, 1, 1);
[root1_20_inc] = inc_dec_fun(x1, 1.01, 1, 1);

d_root0_20_inc = ((root0_20_inc - root0) / root0 * 100);
d_root1_20_inc = ((root1_20_inc - root1) / root1 * 100);

%% b) 

[root0_20_dec] = inc_dec_fun(x0, 0.99, 1, 1);
[root1_20_dec] = inc_dec_fun(x1, 0.99, 1, 1);

d_root0_20_dec = ((root0_20_dec - root0) / root0 * 100);
d_root1_20_dec = ((root1_20_dec - root1) / root1 * 100);

%% c) 

[root1_61_3_inc] = inc_dec_fun(x1, 1, 1.01, 1.01);
[root1_61_3_dec] = inc_dec_fun(x1, 1, 0.99, 0.99);

d_root1_61_3_inc = ((root1_61_3_inc - root1) / root1 * 100);
d_root1_61_3_dec = ((root1_61_3_dec - root1) / root1 * 100);

% Felgräns

%% Utskrifter

roots = [root0; root1];
d_roots_20_inc = [d_root0_20_inc; d_root1_20_inc];
d_roots_20_dec = [d_root0_20_dec; d_root1_20_dec];
d_roots_61_3_dec = [0; d_root1_61_3_dec];
d_roots_61_3_inc = [0; d_root1_61_3_inc];

T1 = table(roots, d_roots_20_inc, d_roots_20_dec, d_roots_61_3_dec, d_roots_61_3_inc, 'VariableNames', {'Rötter', '20 + 1%', '20 - 1%', '61 & 3 + 1%', '61 & 3 - 1%'});
disp('<strong>Tillförlitlighet:</strong>');
fprintf('\n');
disp(T1);

function [root] = inc_dec_fun(x_init, delta_20, delta_61, delta_3)
    f = @(x) 61 * delta_61 * x - ((x^2 + x + 0.03) / (3 * delta_3 * x + 1))^7 - 20 * delta_20 * x * exp(-x);
    options = optimoptions(@fsolve, 'Display', 'off', 'MaxIter', 1e4);
    root = fsolve(f, x_init, options);
end
