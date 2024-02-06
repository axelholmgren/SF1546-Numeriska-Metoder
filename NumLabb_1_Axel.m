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

%% Uppgift 5 
