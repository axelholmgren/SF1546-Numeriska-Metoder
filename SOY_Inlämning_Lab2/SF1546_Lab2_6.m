% SF1546 - VT24 
% Laboration 2 
% Nikolaos Timoudas % Erik Simert

close all 
clear all
clc

%% 6. Ickelinjärt cirkel-problem 

% a) 

punkt1 = [10 , 10];
punkt2 = [12, 2];
punkt3 = [3, 8];

[X_a, Y_a, R_a] = non_lin_newton(punkt1, punkt2, punkt3);

part_a = {'a)'};
T_a = table(part_a, X_a, Y_a, R_a, 'VariableNames', {'Svar', 'X', 'Y', 'R'});

%% b) 

phi = linspace(0, 2*pi, 1e3);
x_a_cirkel = X_a + R_a .* cos(phi); 
y_a_cirkel = Y_a + R_a .* sin(phi);

%% c) 

[X_c, Y_c, R_c] = non_lin_newton_2(punkt1, punkt2, punkt3);

x_c_cirkel = X_c + R_c .* cos(phi); 
y_c_cirkel = Y_c + R_c .* sin(phi);

part_c = {'c)'};
T_c = table(part_c, X_c, Y_c, R_c, 'VariableNames', {'Svar', 'X', 'Y', 'R'});

%% d)

punkt4 = [11, 11];
punkt5 = [2, 9];

[X_d, Y_d, R_d] = newton_gauss(punkt1, punkt2, punkt3, punkt4, punkt5);

x_d_cirkel = X_d + R_d .* cos(phi); 
y_d_cirkel = Y_d + R_d .* sin(phi);

part_d = {'d)'};
T_d = table(part_d, X_d, Y_d, R_d, 'VariableNames', {'Svar', 'X', 'Y', 'R'});

%% e)

[X_e, Y_e, R_e] = newton_gauss_2(punkt1, punkt2, punkt3, punkt4, punkt5);

x_e_cirkel = X_e + R_e .* cos(phi); 
y_e_cirkel = Y_e + R_e .* sin(phi);

figure; hold on; grid on;
p1 = plot(punkt1(1), punkt1(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
p2 = plot(punkt2(1), punkt2(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
p3 = plot(punkt3(1), punkt3(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
p4 = plot(punkt4(1), punkt4(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
p5 = plot(punkt5(1), punkt5(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

cirkel_a = plot(x_a_cirkel, y_a_cirkel, 'b');
mitt_a = plot(X_a, Y_a, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

%cirkel_c = plot(x_c_cirkel, y_c_cirkel);
%mitt_c = plot(X_c, Y_c, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

cirkel_d = plot(x_d_cirkel, y_d_cirkel, 'g');
mitt_d = plot(X_d, Y_d, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'g');

cirkel_e = plot(x_e_cirkel, y_e_cirkel, 'r');
mitt_e = plot(X_e, Y_e, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

title('Lösningar till systemet')
legend([cirkel_a, mitt_a, cirkel_d, mitt_d, cirkel_e, mitt_e], ...
    'Cirkel a', 'Mittpunkt a', 'Cirkel d', 'Mittpunkt d', 'Cirkel e', 'Mittpunkt e');

set(gca,'FontSize',16); set(gca,'FontName','times');
hold off

part_e = {'e)'};
T_e = table(part_e, X_e, Y_e, R_e, 'VariableNames', {'Svar', 'X', 'Y', 'R'});

%% g)

[X_g, Y_g, R_g] = linear_newton(punkt1, punkt2, punkt3);
part_g = {'g)'};
T_g= table(part_g, X_g, Y_g, R_g, 'VariableNames', {'Svar', 'X', 'Y', 'R'});

%% h) 

[X_h, Y_h, R_h] = linear_newton_2(punkt1, punkt2, punkt3);
part_h = {'h)'};
T_h = table(part_h, X_h, Y_h, R_h, 'VariableNames', {'Svar', 'X', 'Y', 'R'});

%% Tabell med resultat  

T_comb = [T_a; T_c; T_d; T_e; T_g; T_h];
disp(T_comb);
 
%% Non - linear Newton 

function [X, Y, R] = non_lin_newton(punkt1, punkt2, punkt3)
    
    x1 = punkt1(1); y1 = punkt1(2);
    x2 = punkt2(1); y2 = punkt2(2);
    x3 = punkt3(1); y3 = punkt3(2);

    % Initiala gissningar 
    X = (x1 + x2 + x3) / 3;
    Y = (y1 + y2 + y3) / 3;
    R = sqrt((x1 - X)^2 + (y1 - Y)^2);

    tol = 1e-6; % Tolerans för konvergens
    n_max = 1e3; % Maximalt antal iterationer

    for i = 1:n_max
        f = [(x1 - X)^2 + (y1 - Y)^2 - R^2;
            (x2 - X)^2 + (y2 - Y)^2 - R^2;
            (x3 - X)^2 + (y3 - Y)^2 - R^2;];

        J = [-2 * (x1 - X), -2 * (y1 - Y), -2 * R;
            -2 * (x2 - X), -2 * (y2 - Y), -2 * R;
            -2 * (x3 - X), -2 * (y3 - Y), -2 * R;];

        t = J \ (-f);
        X = X + t(1);
        Y = Y + t(2);
        R = R + t(3);

        if max(abs(t)) < tol
            break;
        end
    end
end

%% Non - linear Newton for *3

function [X, Y, R] = non_lin_newton_2(punkt1, punkt2, punkt3)
    
    x1 = punkt1(1); y1 = punkt1(2);
    x2 = punkt2(1); y2 = punkt2(2);
    x3 = punkt3(1); y3 = punkt3(2);

    % Initiala gissningar  
    X = (x1 + x2 + x3) / 3;
    Y = (y1 + y2 + y3) / 3;
    R = sqrt((x1 - X)^2 + (y1 - Y)^2);

    tol = 1e-6; % Tolerans för konvergens
    n_max = 1e3; % Maximalt antal iterationer

    for i = 1:n_max
        f = [(x1 - X)^2 + (y1 - Y)^2 - R^2;
            (x2 - X)^2 + (y2 - Y)^2 - R^2;
            3*(x3 - X)^2 + 3*(y3 - Y)^2 - 3*R^2;];

        J = [-2 * (x1 - X), -2 * (y1 - Y), -2 * R;
            -2 * (x2 - X), -2 * (y2 - Y), -2 * R;
            -6 * (x3 - X), -6 * (y3 - Y), -6 * R;];

        t = J \ (-f);
        X = X + t(1);
        Y = Y + t(2);
        R = R + t(3);

        if max(abs(t)) < tol
            break;
        end
    end
end

%% Newton - Gauss function 

function [X, Y, R] = newton_gauss(punkt1, punkt2, punkt3, punkt4, punkt5)
    
    x1 = punkt1(1); y1 = punkt1(2);
    x2 = punkt2(1); y2 = punkt2(2);
    x3 = punkt3(1); y3 = punkt3(2);
    x4 = punkt4(1); y4 = punkt4(2);
    x5 = punkt5(1); y5 = punkt5(2);

    % Initiala gissningar  
    X = (x1 + x2 + x3 + x4 + x5) / 5;
    Y = (y1 + y2 + y3 + y4 + y5) / 5;
    R = sqrt((x1 - X)^2 + (y1 - Y)^2);

    tol = 1e-6; % Tolerans för konvergens
    n_max = 1e3; % Maximalt antal iterationer

    for i = 1:n_max
        f = [(x1 - X)^2 + (y1 - Y)^2 - R^2;
            (x2 - X)^2 + (y2 - Y)^2 - R^2;
            (x3 - X)^2 + (y3 - Y)^2 - R^2;
            (x4 - X)^2 + (y4 - Y)^2 - R^2;
            (x5 - X)^2 + (y5 - Y)^2 - R^2;];
            
        J = [-2 * (x1 - X), -2 * (y1 - Y), -2 * R;
            -2 * (x2 - X), -2 * (y2 - Y), -2 * R;
            -2 * (x3 - X), -2 * (y3 - Y), -2 * R;
            -2 * (x4 - X), -2 * (y4 - Y), -2 * R;
            -2 * (x5 - X), -2 * (y5 - Y), -2 * R;];

        t = J \ (-f);
        X = X + t(1);
        Y = Y + t(2);
        R = R + t(3);

        if max(abs(t)) < tol
            break;
        end
    end
end

%% Newton - Gauss function for *3 in third function 

function [X, Y, R] = newton_gauss_2(punkt1, punkt2, punkt3, punkt4, punkt5)
    
    x1 = punkt1(1); y1 = punkt1(2);
    x2 = punkt2(1); y2 = punkt2(2);
    x3 = punkt3(1); y3 = punkt3(2);
    x4 = punkt4(1); y4 = punkt4(2);
    x5 = punkt5(1); y5 = punkt5(2);

    % Initiala gissningar  
    X = (x1 + x2 + x3 + x4 + x5) / 5;
    Y = (y1 + y2 + y3 + y4 + y5) / 5;
    R = sqrt((x1 - X)^2 + (y1 - Y)^2);

    tol = 1e-6; % Tolerans för konvergens
    n_max = 1e3; % Maximalt antal iterationer

    for i = 1:n_max
        f = [(x1 - X)^2 + (y1 - Y)^2 - R^2;
            (x2 - X)^2 + (y2 - Y)^2 - R^2;
            3*(x3 - X)^2 + 3*(y3 - Y)^2 - 3*R^2;
            (x4 - X)^2 + (y4 - Y)^2 - R^2;
            (x5 - X)^2 + (y5 - Y)^2 - R^2;];
            
        J = [-2 * (x1 - X), -2 * (y1 - Y), -2 * R;
            -2 * (x2 - X), -2 * (y2 - Y), -2 * R;
            -6 * (x3 - X), -6 * (y3 - Y), -6 * R;
            -2 * (x4 - X), -2 * (y4 - Y), -2 * R;
            -2 * (x5 - X), -2 * (y5 - Y), -2 * R;];

        t = J \ (-f);
        X = X + t(1);
        Y = Y + t(2);
        R = R + t(3);

        if max(abs(t)) < tol
            break;
        end
    end
end

%% Linear system 

function [X, Y, R] = linear_newton(punkt1, punkt2, punkt3)
    
    x1 = punkt1(1); y1 = punkt1(2);
    x2 = punkt2(1); y2 = punkt2(2);
    x3 = punkt3(1); y3 = punkt3(2);
    
    b1 = (x1^2 + y1^2);
    b2 = (x2^2 + y2^2);
    b3 = (x3^2 + y3^2);
    b = [b1; b2; b3];
    A = [1, x1, y1;
        1, x2, y2; 
        1, x3, y3];
    c = A\b;
    X = c(2)/2;
    Y = c(3)/2;
    R = sqrt(c(1) + (c(2)^2 + c(3)^2)/4);

end 

%% Linear system med *4 

function [X, Y, R] = linear_newton_2(punkt1, punkt2, punkt3)
    
    x1 = punkt1(1); y1 = punkt1(2);
    x2 = punkt2(1); y2 = punkt2(2);
    x3 = punkt3(1); y3 = punkt3(2);
    
    b1 = (x1^2 + y1^2);
    b2 = (x2^2 + y2^2);
    b3 = (x3^2 + y3^2);
    b = [b1; b2; b3];
    A = [1, x1, y1;
        1, x2, y2; 
        4, 4*x3, 4*y3];
    c = A\b;
    X = c(2)/2;
    Y = c(3)/2;
    R = sqrt(c(1) + (c(2)^2 + c(3)^2)/4);

end 