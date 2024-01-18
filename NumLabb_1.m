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



%% Uppg 3 

function y = f(x)
y = 61*x - ((x^2+x+0.03)/(3*x+1))^7 - 20*x*e^(-x);

end 



