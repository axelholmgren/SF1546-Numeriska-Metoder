% SF1546 - VT24
% Laboration 1 - Uppgift 1
% Nikolaos Timoudas & Erik Simert

close all
clear all
clc

%% U1 Linjärt ekvationssystem
A = [1 2 3 0;
    0 4 5 6;
    1 1 -1 0;
    1 1 1 1];
b = [7 6 5 4]';
x = A\b;

cond = norm(A)*(norm(inv(A)));
r = b - A*x
norm(r)/cond
% decimal64 e_mach = 5*10^-15


