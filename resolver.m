function [MY, MF] = resolver(sis_eqs, t, Y0, metodo)
    % Retorna a solução numérica para um sistema de equações diferenciais

    % Define passo e numero de iteracoes
    h = t(2) - t(1);
    iteracoes = length(t)-1;

    % Inicialização da matriz na qual cada coluna será o vetor Y calculado
    % em um instante de tempo
    MY = zeros(length(Y0), length(t));
    MY(:, 1) = Y0;
    
    % Inicialização da matriz na qual cada coluna será a derivada
    % do vetor Y calculada em um instante de tempo
    MF = zeros(length(Y0), iteracoes);

    % Resolve conforme o método escohido
    
    % Método de Euler
    if metodo == "euler"
        for i=1:iteracoes

            Y = MY(:,i);
            K = sis_eqs(t, Y);
            Y = MY(:,i) + h * K;
            MY(:,i+1) = Y;
            MF(:,i) = K;

        end

    % Método de Runge-Kutta de ordem 2
    elseif metodo == "rk2"
        for i=1:iteracoes

            Y = MY(:,i);
            K1 = sis_eqs(t, Y);
            Y = MY(:,i) + h/2 * K1;
            K2 = sis_eqs(t, Y);
            Y = MY(:,i) + h * K2;
            MY(:,i+1) = Y;
            MF(:,i) = K2;

        end

    % Método de Runge-Kutta de ordem 4
    elseif metodo == "rk4"
        for i = 1:iteracoes

            Y = MY(:,i);
            K1 = sis_eqs(t, Y);
            Y = MY(:,i) + h/2 * K1;
            K2 = sis_eqs(t, Y);
            Y = MY(:,i) + h/2 * K2;
            K3 = sis_eqs(t, Y);
            Y = MY(:,i) + h * K3;
            K4 = sis_eqs(t, Y);
            Y = MY(:,i) + h/6 * (K1 + 2*K2 + 2*K3 + K4);
            MY(:, i+ 1) = Y;
            MF(:, i) = 1/6 * (K1 + 2*K2 + 2*K3 + K4);

        end

    else
        disp("Metodo inválido: " + metodo)    
    end

end
