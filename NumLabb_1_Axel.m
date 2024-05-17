function secant_method()
    % Definiera funktionen
    f = @(x) 61 * x - ((x.^2 + x + 0.03) ./ (3 * x + 1)).^7 - 20 * x .* exp(-x);

    % Initiala uppskattningar
    x0 = 0.0;
    x1 = 0.6;

    % Toleranskrav
    tol = 1e-8;
    
    % Beräkna roten
    [root, iterations, rel_errors] = secant(f, x0, x1, tol);
    
    % Beräkna konvergensordningen
    convergence_order = log(rel_errors(end) / rel_errors(end-1)) / log(rel_errors(end-1) / rel_errors(end-2));
    
    % Visa resultat
    fprintf('Rot: %.12f\n', root);
    fprintf('Antal iterationer: %d\n', iterations);
    fprintf('Konvergensordning: %.3f\n', convergence_order);
end

function [root, iterations, rel_errors] = secant(f, x0, x1, tol)
    % Lista för relativa fel
    rel_errors = [];
    
    iterations = 0;
    
    while true
        % Beräkna ny punkt
        fx0 = f(x0);
        fx1 = f(x1);
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);

        % Relativt fel
        rel_error = abs((x2 - x1) / x2);
        rel_errors = [rel_errors; rel_error];
        
        % Kontrollera tolerans
        if rel_error < tol
            break;
        end

        % Flytta fram punkterna
        x0 = x1;
        x1 = x2;
        iterations = iterations + 1;
    end
    
    root = x1;
end


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

% A) Interpolationspolynom över alla punkter 
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

    
% B) Styckvis linjär interpolation genom samma punkter 
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


% C) Splines-approximation genom alla punkter 
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

    
% D) Andragradspolynom mellan 1 juni- 1 augusti
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


% E) Minsta kvadratanpassat andragradspolynom 1 april- 1 augusti
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


% F) Minstakvadratanpassat andragradspolynom 1 januari- 31 december
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


% G) Minstakvadrat med sinusanpassning
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


% Uppgiftssvar:
% a)
% I vilken ansats behövdes flest koefficienter beräknas (totalt över hela 
% intervallet)?
% Hur många koefficienter är det?
% I kubisk styckvis interpolation är den interpolerade funktionen en serie 
% kubiska polynom. För att hitta koefficienterna var varje polynom behöves
% krävs det n-1 funktioner då vi har samma antal okända. Varje polynom har 
% sedan 4 koefficienter vilket ger 4(n-1) där n är antalet datapunkter 
% (12). Detta ger i detta fall 44 koefficienter

% b)
% De tre anpassningar som bara behöver tre koefficienter är 
% D) Andragradspolynom mellan 1 juni- 1 augusti
% E) Minsta kvadratanpassat andragradspolynom 1 april- 1 augusti
% F) Minstakvadratanpassat andragradspolynom 1 januari- 31 december
% G) Minstakvadrat med sinusanpassning
% då dessa enbart är ett polynom av grad 2 förutom i G) då denna ger 
% koefficienter amplitud, fasförskjutning och skiftnign i y-led.

% c) 
% För att beräkna den längsta dagen ges att alternativ 
% D) Andragradspolynom mellan 1 juni- 1 augusti. Detta för att det ger en 
% bra anpassning för intervallet där den längsta dagen infaller samt att
% den ger effektiv beräknignsmetod där det enbart krävs 3 koefficienter.

% d)
% Julafton infaller vid dag 358 av 365 vilket sätter krav på anpassningen
% vid ändpunkterna. Här ger
% B) Styckvis linjär interpolation genom samma punkter 
% C) Splines-approximation genom alla punkter 
% G) Minstakvadrat med sinusanpassning
% Här har B) inte nog med prestanda då den enbart är en styckvis 
% interpolation mellan punkterna, den är snabb men inte en bra
% anpassning. Här ger spline den bästa prestandan men har 44 koefficienter
% medan G) ger relativt bra prestanda men med enbart 3 koefficienter.
% Väljer G) som den bästa anpassningen.

% e)
% En potentiell förbättring är att använda D) Andragradspolynom mellan 
% 1 juni- 1 augusti med modifikationen att använda data mellan 1 nov- 31
% dec där julafton ligger i intervallet

subplot(3,3,8);
    delX = x(11:13);
    delY = y(11:13);
    % Grad 2
    k4 = polyfit(delX,delY,2);
    p2 = polyval(k4,d);
    plot(d,p2);
    hold on;
    plot(x,y, 'b.');
    title("2a gradspolynom (1 nov-31 dec )");
    xlabel("Dagar");
    ylabel("Minuter");
    axis([0,365,0,1300]);
    hold on;

% f) Olika metoder ger olika fördelar men för att täcka in flest
% datapunkter med bra prestanda ges av C) Splines-approximation genom 
% alla punkter. Däremot är den inte lika effektiv som de andra metoderna.
% Det går inte att säga att en är bättre än den andra, det beror på vad som
% pritoriteras.

%% Uppgift 5 Tillförlitlighet

%% Från Uppgift 2:
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
%% 

% Initierar ett figurfönster för att presentera plottar.
figure();

% Definierar funktionens definitionsmängd genom att definitera: x mellan -5
% och 7.5 med steglängd 1e-5
x= -5:1e-5:7.5;

% Definierar funktionen f(x)
f = @(x) 61*x - ((x.^2 + x + 0.03)./(3*1.01*x + 1)).^7 - 20*x.*exp(-x);


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

disp("  _____Uppg 5: _________ ")

% Derivera funktionen f över x
df_dx= @(x) 61 - 7*((x.^2+x+0.03)./(3*1.01*x+1)).^6*( ( (3*1.01*x+1)*(2*x+1) - (x.^2+x+0.03)*(3*1.01) )./( (3*1.01*x+1)^2 ) ) - 20*(exp(-x)-x.*exp(-x));

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


% a) 
% Den minsta postitiva roten blir 2,614777619×10^-15 mindre vilket motsvarar
% en procentuell ändring på ca 1 %.
% Den största postitiva roten blir 6,4899*10^-6 större vilket motsvarar en
% procentuell ökning på ca 1%. 
% Rötterna påverkas ungefär lika mycket procentuellt sett men påverkas åt
% olika håll. Varför ändringen blir som den blir beror på att kosntanten
% står framför en negativ exponent vilket betyder att den har exponentiellt
% mindre påverkan för högre värden på x

% b)
% Den minsta postivia roten blir 2,5893914285×10^-15 större vilket motsvara
% en procentuell ändring på ca 1% 
% Den största postivita roten blir 6,4898*10^-6 mindre vilket motsvara en
% procentuell ändring på ca 1 %. Det går här att se att rötterna ändra
% ungefär lika mycket som i a) fast med omvänt tecken. Detta beror på samma
% förklaring som i a)

% c)
% Felfortplantningen ger att Ef= 0.1973846874 genom den allmänna
% felfortplantningsformeln och ger felgränserna: 3.397062994660476 (+-)
% 0.1973846874

%% 7. Integral med förbehandling





