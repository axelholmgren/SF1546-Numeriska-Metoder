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
 

%b)
% Initierar ett figurfönster för att presentera plottar.
figure();

% Definierar funktionens definitionsmängd genom att definitera: x mellan -5
% och 7.5 med steglängd 1e-5
x= -5:1e-5:7.5;

% Definierar funktionen f(x)
f = @(x) 61*x - ((x.^2 + x + 0.03)./(3*x + 1)).^7 - 20*x.*exp(-x);


% Ger funktionens värdemängd, y är resultatet av x genom funktionen.
% Skapar en vektor av samma längd
y = f(x);

% Skapar 5 subplottar, 4 för rötter och en för funktionen i helhet.
subplot(3, 2, 1);
    plot(x, y);
        title('Rot 1');
        xlim([-1.15, -1.1]);
        ylim([-1, 1]);
        yline(0, 'color' ,'k', 'linestyle', ':', 'LineWidth', 1);

subplot(3, 2, 2);
    plot(x, y);
        title('Rot 2');
        xlim([-0.30, -0.28]);
        ylim([-1, 1]);
        yline(0, 'color' ,'k', 'linestyle', ':', 'LineWidth', 1);

subplot(3, 2, 3);
    plot(x, y);
        title('Rot 3');
        xlim([-0.001, 0.001]);
        ylim([-0.05, 0.05]);
        yline(0, 'color' ,'k', 'linestyle', ':', 'LineWidth', 1);

subplot(3, 2, 4);
    plot(x, y);
        title('Rot 4');
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


% Från plottarna can rötterna grafiskt lokaliseras i följande:
x1=-1.15;
x2=-0.29;
x3=1e-2;
x4=6.34;


% c) Ta fram en av rötterna med Newtons Metod.
disp("  _____Uppg 2 c): _________ ")

% Derivera funktionen f över x
df_dx= @(x) 61 - 7*((x.^2+x+0.03)./(3*x+1)).^6*( ( (3*x+1)*(2*x+1) - (x.^2+x+0.03)*(3) )./( (3*x+1)^2 ) ) - 20*(exp(-x)-x.*exp(-x));

%Initierar variabler x=input av funktion
% Val av x1,x2,x3,x4
x=x4;
t=1;
x_forra=1;
rel_fel=1;
format short 
disp("    ___________________Newton Raphson_____________________")
disp("      x       fx       dfx_dx       t       lin       kvad");
%Iterativ metod med Newton raphson
while abs(rel_fel)>1e-8
    fx=f(x);
    dfx_dx=df_dx(x);

    %Newton Raphsons metod
    t_forra=t;
    t=fx/dfx_dx;
    
    lin=t/t_forra;
    kvad=t/(t_forra^2);
    
    disp([x fx dfx_dx t lin kvad]);
    x_forra=x;
    x=x-t;
    rel_fel=(x-x_forra)/x;
end;
format long e
rot_nr=x

rot_1_nr=-1.115141590525118e+00;
rot_2_nr=-2.910130809545486e-01;
rot_3_nr=5.334146342066178e-13;
rot_4_nr=6.397062994660476e+00;
rotter=[rot_1_nr rot_2_nr rot_3_nr rot_4_nr];

format short 

disp(["Rot 1" "Rot 2" "Rot 3" "Rot 4"])
disp([ rotter ])




% d)
% Kvadratisk konvergens definieras som hastigheten med vilken felet minskar
% det vill säga att en lösningsmetod har kvadratisk konvergens om felet
% minskar kvadratiskt, dvs att felet är proportionellt med kvadraten av det
% föregående felet. Alltså är det mått på hur snabbt metoden närmar sig
% rätt svar.
% 
% lim (n → ∞) |e_(n+1)| / |e_n|^2 = Konstant


%% e)
disp("    ____Uppg 2 e) :_____  ")

%Upprepar Newton raphsons och plottar differensen: linjärt och kvadratisk 
% konvergens för att se om det finns något sådant samband


%Initierar variabler x=x4 för största roten
x=x4;
t=1;
x_forra=5;
rel_fel=1;
t=1;

disp("    ___________Newton Raphson_____________")
disp("      x       fx       lin       kvad");

%Upprepar iterativ metod med Newton raphson
%För att läsa ut konvergens-konstanten
while abs(rel_fel)>1e-8
    fx=f(x);
    dfx_dx=df_dx(x);

    %Newton Raphsons metod
    t_forra=t;
    t=fx/dfx_dx;
    
    lin=t/t_forra;
    kvad=t/(t_forra^2);
    
    disp([x fx lin kvad]);
    x_forra=x;
    x=x-t;
    rel_fel=(x-x_forra)/x;
end;

rot_NR=x

%Vi kan se att vi har kvadratisk konvergens då t/t_förra^2 är så gott som
% konstant

%% Uppg 3 
disp("      _______ Uppg 3 a): _______")

% Funktionen f definieras
f = @(x) 61*x - ((x.^2 + x + 0.03)./(3 * x + 1 )).^7 - 20*x.*exp(-x);



disp("   ______Sekantmetoden_______");

%Initialvärden
x0=0;
x1=1;

x=x1;
x_forra=x0; %x0
f_forra=f(0); %f0
t_forra=1;
t_foforra=1;
disp("   x   f   korrektionsterm konvvergens") 
while abs(x-x_forra)>1e-8;

    f_n=f(x); %Beräknar nuvarande funktionsvärdet

    t_n=f_n*(x-x_forra)/(f_n-f_forra);


    disp([x_forra f_forra])


    disp([x f_n t_n]) 
    x_forra=x;

    k=t_n/(t_forra*t_foforra);

    f_forra=f_n; 
    t_forra=t_n
    t_foforra=t_forra;
    x=x-t_n; 
    

end; 
rot=x

% b)

% Konvergensordning ca 1.6 (dvs mellan linjär och kvadratisk)


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
    title("2a gradspolynom (1 jun-1 aug)");
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
    title("2a gradspolynom (1 apr-1 sept)");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0,365,0,1300]);
    hold on;


% F) 
subplot(3,3,6);
    % Grad 2
    k6 = polyfit(x,y,2);
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
    w = 2*pi/365;
    c1 = ones(13,1);
    c2 = cos(w*x);
    c3 = sin(w*x);
    C = [c1,c2,c3];

    %beräknar minsta kvadrat c

    c = (C'*C)^(-1)*(C'*y);

    p5 = c(1)+c(2)*cos(w*d)+c(3)*sin(w*d);

    plot(d,p5);
    hold on;
    plot(x,y,'b.');
    title("Minstakvadratanpassad med sinus");
    xlabel("Dagar")
    ylabel("Minuter");
    axis([0,365,0,1300]);


    
    





