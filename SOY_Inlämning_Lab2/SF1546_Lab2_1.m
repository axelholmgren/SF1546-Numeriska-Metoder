% SF1546 - VT24 
% Laboration 2 
% Nikolaos Timoudas % Erik Simert

clear 
close all
clc 

%% 1. Numerisk integration: Eulers metod 

% a) 

%Differentialekvation 
f = @(x,y) -(1/6+(pi.*sin(pi.*x))./(1.6-cos(pi.*x))).*y;

%Steglängder
h = 0.5;
steg_h = [h, h/2, h/4, h/8, h/16, h/32, h/64, h/128];

%Konstanter 
x0 = 0;
x1 = 4;

%Eulers metod 
for i = 1:length(steg_h)
    n = (x1-x0)./steg_h(i); %antal steg 
    xVec = (x0:steg_h(i):x1); %vektor för x
    y = zeros(1, n+1); %allokera minne till y 
    y(1) = 2.5; %begynnelsevärde 

    for j=1:n
        y(1+j)=y(j)+steg_h(i)*f(xVec(j),y(j)); %Euler
    end
    y_4(i) = y(end); %lagra y(4) för de olika steglängderna  
    
    if i > 1
        error(i) = y_4(i) - y_4(i-1);
    end
    % Plot 
    plot(xVec, y, 'DisplayName', ['h = ', num2str(steg_h(i))]);
    hold on
end

legend('Location', 'northeast');
xlabel('x');
ylabel('y(x)');
title('Eulers Metod');
grid on;
set(gca,'FontSize',16);
set(gca,'FontName','times');
hold off 


T1 = table(steg_h', y_4', error', 'VariableNames', {'Steglängd', 'y(4)', 'Felgräns'});
disp(T1)


%% Felgräns 

% b) 

% Säker decimal då h / 128 

%error = abs(y_4(end) - y_4(end-1))


