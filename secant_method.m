function sekantmetoden()
    % Definiera funktionen enligt det givna uttrycket
    f = @(x) 61 * x - ((x.^2 + x + 0.03) ./ (3 * x + 1)).^7 - 20 * x .* exp(-x);

    % Skapa en array av x-värden för att plotta
    x = linspace(-5, 5, 1000);  % Justera intervallet efter behov
    y = arrayfun(f, x);  % Utvärdera funktionsvärden

    % Plotta funktionens beteende
    figure;  % Skapa en ny figur
    plot(x, y, 'b-', 'LineWidth', 1.5);
    xlabel('x');
    ylabel('f(x)');
    title('Plott av funktionen f(x)');
    grid on;

    % Valfri: Rita en horisontell linje vid y=0 för referens
    hold on;
    yline(0, 'r--');
    hold off;

    % Välj initiala gissningar för sekantmetoden
    x0 = -1.14;  % Justera till ett intervall nära förväntade rötter
    x1 = -1.15;  % Justera på liknande sätt

    % Konvergenstolerans
    tol = 1e-8;

    % Hitta roten med sekantmetoden
    [rot, iterationer, rel_fel] = sekant(f, x0, x1, tol);

    % Visa resultaten
    fprintf('Rot: %.12f\n', rot);
    fprintf('Iterationer: %d\n', iterationer);

    % Beräkning av konvergensordning
    if length(rel_fel) >= 3
        konvergens_ordning = log(rel_fel(end) / rel_fel(end-1)) / log(rel_fel(end-1) / rel_fel(end-2));
        fprintf('Konvergensordning: %.3f\n', konvergens_ordning);
    else
        fprintf('Inte tillräckligt med data för att beräkna konvergensordning.\n');
    end
end

function [rot, iterationer, rel_fel] = sekant(f, x0, x1, tol)
    % Lista för relativa fel
    rel_fel = [];
    
    iterationer = 0;
    
    while true
        % Utvärdera funktionsvärden
        fx0 = f(x0);
        fx1 = f(x1);

        % Undvik division med noll
        if fx1 - fx0 == 0
            fprintf('Sekantmetoden misslyckades: Division med noll.\n');
            return;
        end

        % Beräkna ny punkt
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);

        % Relativt fel
        rel_fel = [rel_fel; abs((x2 - x1) / x2)];
        
        % Kontrollera tolerans
        if rel_fel(end) < tol
            break;
        end

        % Skjut fram punkterna
        x0 = x1;
        x1 = x2;
        iterationer = iterationer + 1;
    end
    
    rot = x1;
end
