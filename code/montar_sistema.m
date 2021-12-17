function sis_eqs = montar_sistema(proporcao_m, tracao)
    % Monta sistema de equações diferenciais conforme o caso escolhido.
    % Retorna function handle com t e o vetor Y como argumentos.

    % Constantes
    mtotal = 1939;
    L = 2.95;
    I = 1;
    w = 10;
    mi = 0.42;
    beta = 0.02;
    g = 9.8;

    % Distribuição de massas
    m = proporcao_m*mtotal;
    m1 = (1 - proporcao_m)*mtotal;

    % Forças atuantes nos pneus
    if tracao == "traseira"
        F = beta*m*g;
        F1 = mi*m1*g;
    elseif tracao == "dianteira"
        F = -mi*m*g;
        F1 = -beta*m1*g;
    elseif tracao == "quatro rodas"
        F = -mi*m*g;
        F1 = mi*m1*g;
    else
        disp("Tração inválda: " + tracao)
    end
    
    % Monta o sistema de equações
    sis_eqs = @(t, Y) [...
        Y(2);
        -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
        + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
        Y(4);
        -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
        - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
        /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
    ];      
    
end