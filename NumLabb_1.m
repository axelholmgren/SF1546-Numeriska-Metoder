%______________Numme Labb 1_________________
% Uppg 1 - Linjärt ekvationssystem

% 1a)

%Vi ska lösa ekvationssystemet Ax=b
%Deffinierar matris A
A=[1 2 3 0;
  0 4 5 6;
  1 1 -1 0;
  1 1 1 1]

b=[7 6 5 4]'

%Definierar x-vektorn
%x=[x1 x2 x3 x4]'

x=A\b

% 1b)
%Beräknar residualvektorn r
r=b-A*x

% 1c)
%residualvektorn blir inte exakt 0 pga att 
%vektordivisionen inte beräknar värdena exakt utan 
%dessa approximeras, och har inte med fullständiga 
%värdet och därmed blir inte b exakt lika med A*x

%% Uppg 2
%% 
%b)
%Initiate a figure-window for presenting our plots
figure();

%Define the domain of the function (definitionsmängd) by defining the scope
%of x ; between -5 and +7.5 with steps of size order 1e-5
x= -5:1e-5:7.5;

%Defines the function f(x)
f = @(x) 61*x - ((x.^2 + x + 0.03)./(3*x + 1)).^7 - 20*x.*exp(-x);

%This gives the range of the function, y is the output when x is run
%through our function - creates a vector of equal length
y = f(x);

%Creates 5 subplots, 4 for all the roots, and one for the full function
subplot(3, 2, 1);
    plot(x, y);
        title('1:st root');
        xlim([-1.15, -1.1]);
        ylim([-1, 1]);
        yline(0, 'color' ,'k', 'linestyle', ':', 'LineWidth', 1);

subplot(3, 2, 2);
    plot(x, y);
        title('2:nd root');
        xlim([-0.30, -0.28]);
        ylim([-1, 1]);
        yline(0, 'color' ,'k', 'linestyle', ':', 'LineWidth', 1);

subplot(3, 2, 3);
    plot(x, y);
        title('3:d root');
        xlim([-0.001, 0.001]);
        ylim([-0.05, 0.05]);
        yline(0, 'color' ,'k', 'linestyle', ':', 'LineWidth', 1);

subplot(3, 2, 4);
    plot(x, y);
        title('4:th root');
        xlim([6.35, 6.45]);
        ylim([-2, 2]);
        yline(0, 'color' ,'k', 'linestyle', ':', 'LineWidth', 1);

subplot(3, 2, 5:6);
    plot(x, y);
        title('y=f(x)');
        xlim([-3, 7]);
        ylim([-300, 300]);        
        yline(0, 'color' ,'k', 'linestyle', ':', 'LineWidth', 1);


%exportgraphics(gcf,'arrayAntennaMeshMew.png','Resolution',300);

%In the graphs we can read out that the roots roughly are located at 
x1=-1.15;
x2=-0.29;
x3=1e-2;
x4=6.34;


% c) Decide one of the roots accuratly with Newtons Method

%Derative of the function f over x
df_dx= @(x) 61 - 7*((x.^2+x+0.03)./(3*x+1)).^6*( ( (3*x+1)*(2*x+1) - (x.^2+x+0.03)*(3) )./( (3*x+1)^2 ) ) - 20*(exp(-x)-x.*exp(-x));

%Initierar variabler x=input of function
x=x1
t=1
x_prev=1
rel_err=1

%Iterative method using Newton raphson
while abs(rel_err)>1e-8
    fx=f(x);
    dfx_dx=df_dx(x);

    %Newton Raphsons method
    g=t;
    t=fx/dfx_dx;
    
    lin=t/g;
    quad=t/(g^2);
    
    disp([x fx dfx_dx t lin quad]);
    x_prev=x;
    x=x-t;
    rel_err=(x-x_prev)/x;
end;

root_NR=x



%% Uppg 3 

% Funktion
f = @(x) 61*x - ((x.^2 + x + 0.03)./(3 * x + 1 )).^7 - 20*x.*exp(-x);



disp("Sekantmetoden");
x0=0;
x1=1;

x_n=x1
x_n_prev=x0;

f_n_prev=f(0) %f0

while abs(x_n-x_n_prev)>5e-5;

    f_n=f(x_n)

    t_n=f_n*(x_n-x_n_prev)/(f_n-f_n_prev)
    disp("x f korr konv") 
    disp([x_n_prev f_n_prev])


    disp([x_n f_n t_n]) 
    x_n_prev=x_n;

    f_n_prev=f_n; 
    x_n=x_n-t_n; 

end; 
rot=x_n



