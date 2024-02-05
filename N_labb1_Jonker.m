%Laboration 1 - Jonathan Fransson
%% -----Uppgift 1-----
%a)

A = [1 2 3 0; 0 4 5 6; 1 1 -1 0; 1 1 1 1];
b = [7 6 5 4]';

x = A\b
%% 
%b)

A = [1 2 3 0; 0 4 5 6; 1 1 -1 0; 1 1 1 1];
b = [7 6 5 4]';

r = b - A*x

%% 
%c)
%På grund av att alla komponenterna i 'x' i detta fallet inte består av
%rationella tal kommer Matlab behöva göra avrundningar. Dessa avrundningar
%är av rangordningen e^-16 eftersom de är de minsta/största tal datorn kan
%hantera. Därav erhålles restvärden när vi tar 'r=b-A*x'.

%% -----Uppgift 2-----
%a) Se bifogat dokument

%% 
%b)
figure();

%Variable
x = -3:1e-4:6.5;

%Function
f = @(x) 61*x - ((x.^2 + x + 0.03)./(3*x + 1)).^7 - 20*x.*exp(-x);
func = f(x);

%Plotting
subplot(3, 3, 1);
    plot(x, func);
        ylim([-1, 1]);
        xlim([-1.22, -1.1]);
        title('Root nr1');
        line([-10, 10], [0, 0], 'color' ,'r', 'linestyle', '--');

subplot(3, 3, 2);
    plot(x, func);
        ylim([-0.1, 0.1]);
        xlim([-0.30, -0.288]);
        title('Root nr2');
        line([-10, 10], [0, 0], 'color' ,'r', 'linestyle', '--');

subplot(3, 3, 3);
    plot(x, func);
        ylim([-0.1, 0.1]);
        xlim([-0.001, 0.001]);
        title('Root nr3');
        line([-10, 10], [0, 0], 'color' ,'r', 'linestyle', '--');

subplot(3, 3, 4);
    plot(x, func);
        ylim([-1, 1]);
        xlim([6.35, 6.45]);
        title('Root nr4');
        line([-10, 10], [0, 0], 'color' ,'r', 'linestyle', '--');

subplot(3, 3, 5);
    plot(x, func);
        ylim([-100, 300]);
        xlim([-2, 7]);
        title('Function view');
        line([-10, 10], [0, 0], 'color' ,'r', 'linestyle', '--');


%exportgraphics(gcf,'arrayAntennaMeshMew.png','Resolution',300);

%% 
%c) Newtons metod

%Variables
x1= -1.3;
x2=-0.29;
x3=0.01;
x4=6;
t=1;

x0=1;
err=1;

%Input
%x=x1;
%x=x2;
%x=x3;
x=x4;

%Formatering
format long e
disp('x f(x) fprim(x) korr kvad linj');

%Algoritm
while abs(err)>1e-8
    %Funktion och derivata
    f = 61*x - ((x^2 + x + 0.03)/(3*x + 1))^7 - 20*x*exp(-x);
    fp = 21*(((x^2 + x + 0.03)^7)/((1 + 3*x)^8)) - 7*(2*x + 1)*(x^2 + x + 0.03)^6/((1 + 3*x)^7) - 20*exp(-x) + 20*exp(-x)*x + 61;
    %f = 5*x^2-10*x-2;
    %fp =10*x-10;

    %Algoritm (Newton-Raphson)
    g=t;
    t=f/fp;
    quad=t/g^2; 
    lin=t/g;
    disp([x f fp t quad lin]);
    x0=x;
    x=x-t;
    err=(x-x0)/x;
end;


%Svar
root=x

%% -----Uppgift 3-----
%a) Sekantmetoden

format long e
disp('Sekantmetoden');

% root = -1.115141590525118e+00;  %stämmer
% x0=-1.5;
% x1=-1;

% root = -2.910130809545486e-01;  %stämmer
% x0=-0.290;
% x1=-0.295;

% root = 5.334146342066178e-13; %stämmer
% x0=-0.1;
% x1=0.1;

% root = 6.397062994660478e+00; %stämmer
% x0=5
% x1=7

g0=1;
g1=1;
f0=61*x0-((x0^2 + x0 + 0.03)/(3*x0 + 1))^7-20*x0*exp(-x0);


disp('x f korr konv');
disp([x0 f0]);
while abs((x1-x0)/x1)>1e-8
    f1=61*x1-((x1^2 + x1 + 0.03)/(3*x1 + 1))^7-20*x1*exp(-x1);
    t=f1*(x1-x0)/(f1-f0); %t = 1/fp
    k=t/(g1*g0);
    disp([x1 f1 t k])
    
    x0=x1;
    f0=f1;
    x1=x1-t;
    g0=g1;
    g1=t;
end;

rot=x1
abs_err = abs((rot - root))
%%
%b) Se dokument

%c) Se dokument
 
%d)

format long e
disp('Sekantmetoden');

% root = -1.182695404466818e+00  Stämmer
x0=-1.5;
x1=-1;

% root = -2.935747112778986e-01  Stämmer
% x0=-0.290;
% x1=-0.295;

% root ≈ 3.810232865426745e-12 Stämmer inte
% x0=-0.1;
% x1=0.1;

% root = 6.414684289741624e+00 stämmer inte
% x0=5
% x1=7

xt1=0;
xt2=0;

g0=1;
g1=1;
f0=62*x0-((x0^2 + x0 + 0.04)/(3*x0 + 1))^7-19*x0*exp(-x0);


disp('x f korr konv');
disp([x0 f0]);
while abs((x1-x0)/x1)>1e-8
    f1=62*x1-((x1^2 + x1 + 0.04)/(3*x1 + 1))^7-19*x1*exp(-x1);
    t=f1*(x1-x0)/(f1-f0); %t = 1/fp
    k=t/(g1*g0);
    disp([x1 f1 t k])
    
    xt2=xt1;
    xt1=x0;
    x0=x1;
    f0=f1;
    x1=x1-t;
    g0=g1;
    g1=t;
end;

uppe=log(abs(x1-x0)/abs(x0-xt1));
nere=log(abs(x0-xt1)/abs(xt1-xt2));

q = uppe/nere

rot=x1

%% -----Uppgift 4-----
%a) Se dokument

%b)
format long e

error=1;
x0=0.1;
xt0=1;
xt1=1;
xt2=1;

while abs(error)>1e-4
    f=62*x0-((x0^2 + x0 + 0.04)/(3*x0 + 1))^7-19*x0*exp(-x0);
    xt2=xt1;
    xt1=xt0;
    xt0=x0;
    x1 = (62*x0 - f)/62;
    error = (x1-x0)/x1;
    x0=x1;
end;

uppe=log(abs(x1-xt0)/abs(xt0-xt1));
nere=log(abs(xt0-xt1)/abs(xt1-xt2));

q = uppe/nere

rot = x1

%% -----Uppgift 5-----

A = [1 32 60 91 121 152 182 213 244 274 305 335 365; 375 486 632 795 955 1083 1104 998 844 684 526 396 374]';

x = A(:,1);
y = A(:,2);

days = 1:1:365;

figure();
%A)
subplot(3,3,1);
    k1 = polyfit(x, y, 12); %Runges fenomen då vi går över grad 10
    f1 = polyval(k1, days);
    plot(days, f1);
    hold on;
    plot(x, y, 'r.');
    title("Interpolations polynom genom alla punkter");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0, 365, 0, 1440]);
    hold on;
   
%B)
subplot(3,3,2);
    k2 = interp1(x, y, days, "linear");
    plot(days, k2);
    hold on;
    plot(x, y, 'r.');
    title("Styckvis interpolation genom alla punkter");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0, 365, 0, 1440]);
    hold on;

%C)
%Splines kan inte gissa framtiden. "Tredjegrads interp" innanför datamängden,
%hanterar inte approx utanför mätvärden

subplot(3,3,3);
    k3 = spline(x, y, days);
    plot(days, k3);
    hold on;
    plot(x, y, 'r.');
    title("Splines-approx genom alla punkter");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0, 365, 0, 1440]);
    hold on;

%D)
subplot(3,3,4);
    %Träffar dom sökta punkterna bra, inte resten, för låg  grad
    slicedX = x(6:8);
    slicedY = y(6:8);

    k4 = polyfit(slicedX, slicedY, 2); %2:a grads polynom
    f4 = polyval(k4, days);
    plot(days, f4);
    hold on;
    plot(x, y, 'r.');
    title("2:a grad int. polynom genom juni-aug");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0, 365, 0, 1440]);
    hold on;

%E)
 subplot(3,3,5);
    slicedX2 = x(4:9);
    slicedY2 = y(4:9);
    
    k5 = polyfit(slicedX2, sliced 2Y2, 2); %2:a grads polynom
    f5 = polyval(k5, days);
    plot(days, f5);
    hold on;
    plot(x, y, 'r.');
    title("2:a grad int. polynom genom april-sept");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0, 365, 0, 1440]);
    hold on;

%F)
subplot(3,3,6);
    k6 = polyfit(x, y, 2); %2:a grads polynom
    f6 = polyval(k6, days);
    plot(days, f6);
    hold on;
    plot(x, y, 'r.');
    title("2:a grad int. polynom genom all data");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0, 365, 0, 1440]);
    hold on;

%G)
subplot(3,3,7);
    w = 2*pi/365;
    cv1 = ones(13,1);
    cv2 = cos(w*x);
    cv3 = sin(w*x);
    cMatrix = [cv1 cv2 cv3];
    
    alph = cMatrix'*cMatrix;
    beta = cMatrix'*y;
    c = alph^(-1)*beta;

    f7 = c(1) + c(2)*cos(w*days) + c(3)*sin(w*days);

    plot(days, f7);
    hold on;
    plot(x, y, 'r.');
    title("Sinusanpassning");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0, 365, 0, 1440]);

%a) 
%Den ansats som behöver mest koefficienter är spline funktionen.
%Antalet koefficienter som behövs ges av ekvationen:
%(k+1)n - 2(k-1)
%där 'k' är graden på polynomet och n är antalet datapunkter.
%Graden på polynomet för Matlabs inbyggda 'spline' funktion är 2, antalet
%datapunkter är 365.
%Ansatsen behöver i detta fall alltså 1093 koefficienter

%b)
%Dom ansatser som enbart behövde 3 koefficienter är:
    %Sinusanpassning
    %2:a grad int. polynom genom juni-aug
    %2:a grad int. polynom genom april-sept
    %2:a grad int. polynom genom all data

%c)


%% -----Uppgift 6-----

%Variables
x1= -1.3;
x2=-0.29;
x3=0.01;
x4=6;
t=1;

x0=1;
err=1;

%Input
x=x1;
%x=x2;
%x=x3;
%x=x4;

k19 = 1.03;
k62 = 1;

%Formatering
format long e
%disp('x f(x) fprim(x) korr kvad linj');

%Algoritm
while abs(err)>1e-8
    %Funktion och derivata
    f = (k62*62)*x - ((x^2 + x + 0.04)/(3*x + 1))^7 - (k19*19)*x*exp(-x);
    fp = 21*(((x^2 + x + 0.04)^7)/((1 + 3*x)^8)) - 7*(2*x + 1)*(x^2 + x + 0.04)^6/((1 + 3*x)^7) - (k19*19)*exp(-x) + (k19*19)*exp(-x)*x + (k62*62);
    %f = 5*x^2-10*x-2;
    %fp =10*x-10;

    %Algoritm (Newton-Raphson)
    g=t;
    t=f/fp;
    quad=t/g^2; 
    lin=t/g;
    %disp([x f fp t quad lin]);
    x0=x;
    x=x-t;
    err=(x-x0)/x;
end;

%errBig62=x
%errBig19=x
%errSmall62=x
%errSmall19=x
%bigroot97=x
%bigroot1=x
%smallroot97=x
%smallroot1=x

smallchange = abs(1-(smallroot1/smallroot97))*10^2
bigchange = abs(1-(bigroot97/bigroot1))*10^2

bigRootError = (abs(1-(errBig62/bigroot1)) + abs(1-(errBig19/bigroot1)))*10^2
smallRootError = (abs(1-(errSmall62/smallroot1)) + abs(1-(errSmall19/smallroot1)))*10^2

%a)
%Den minsta roten minskar och den största roten ökar. Den minsta roten minskar med ca 2.5% medans den
%största roten ökar med ca 0.0003%. Detta är pågrund av att termen med
%konstanten 19 framför är en negativ exponent och har därför exponentiellt mindre
%inverkan vid högre värden på x.


%c)
%För den stora roten så blir felet ca 0.55% och för den lilla roten blir
%felet ca 5%. Anledningen för varför det blir större diff på den minsta roten
%än för den största är samma motivering som i a)


%% -----Uppgift 7-----


% x_min = 0;        % Minimum x
% x_max = 1e-4;       % Maximum x
% %x_max = 1.5;
% N = 10000;           % Discretization of x
% x = linspace(x_min, x_max, N);      % x-vector
% 
% % Create integrand
% integrand = (1 - exp(- ((x./3).^3)) ) .* 1./(5*x.^3);
% integrand2 = (-expm1(- ((x./3).^3)) ) .* 1./(5*x.^3);
% integrand3 =(1./(5*x.^3));
% 
% % -------------------------------------------------------------------------
% % 1 - exp(a) for small a is computationally problematic
% % We have that 
% 
% % -------------------------------------------------------------------------
% % expm1(a) = exp(a)-1, with higher precision for small a!
% %
% % 
% figure(1)
% plot(x,integrand,'Linestyle','--','LineWidth',0.9)
% hold on
% %plot(x,integrand3,'LineWidth',0.9)
% %hold on
% plot(x,integrand2,'LineWidth',0.9)
% xlim([x_min x_max])
% %ylim([0 1]);


%Trapets intergration

h = 0.001;
tolerance = 10e-11;
totArea = 0;
refArea = 1;
tailArea = 1.000000000000000e-11;
b=100000;
Etrunk = 1;

while Etrunk>tolerance
    x0 = 0;
    totArea = 0;
    while x0<b
        if x0 < 1e-4
            totArea = totArea + (h/2)*(((1/135)*exp((-x0^3)*1/27))+((1/135)*exp((-(x0+h)^3)*1/27)));
            x0 = x0+h;
            continue;
        end;
        totArea = totArea + (h/2)*((1-exp(-((x0/3)^3)))*1/(5*x0^3)+(1-exp(-(((x0+h)/3)^3)))*1/(5*(x0+h)^3));
        x0 = x0+h;
    end;
    Etrunk = abs(totArea - refArea);
    refArea = totArea;
    h=h/2;
end;



totArea = totArea + tailArea;

totArea

%% -----Uppgift 8-----
format long e
%N = 100000000;           % Discretization of x
x_min = 0.135;        % Minimum x
x_max = 0.138;       % Maximum x
%x = linspace(x_min, x_max, N);      % x-vector

%integrand1 = 159.*exp(-(((23.*x - pi)./0.002).^2));

%figure(1);
%plot(x, integrand1, 'LineWidth', 0.9);


%Quad 
tolerance = 1e-10;
q = quad(@(x) 159*exp(-(((23.*x - pi)./0.002).^2)), x_min, x_max, tolerance)

%Integral
i = integral(@(x) 159*exp(-(((23.*x - pi)./0.002).^2)), x_min, x_max, "AbsTol", tolerance)