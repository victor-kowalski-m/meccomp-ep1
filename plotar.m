function plotar(metodo, t, h, MY, MF, proporcao_m, tracao)

    % Titulo geral do gráfico e da imagem salva
    titulo = tracao + " m " + proporcao_m + " m1 " + (1 - proporcao_m )+ ...
        " " + metodo + " (h = " + h + ")";

    % Plot eixo x
    figx = figure('visible','off');
    yyaxis left
    plot(t, MY(2, :), '-')
    ylabel('Velocidade (m/s)')
    yyaxis right
    plot(t(2:end), MF(2, :), '-')
    ylabel('Aceleração (m/s^2)')
    title("x " + titulo)
    grid()

    % Plot angulo theta
    figtheta = figure('visible','off');
    yyaxis left
    plot(t, MY(3, :), '-')
    ylabel('Posição Angular (rad)')
    yyaxis right
    plot(t, MY(4, :), '-')
    ylabel('Velocidade Angular (rad/s)')
    title("theta " +  titulo)
    grid()

    % Escolhe qual pasta salvar, dependendo da parte
    if proporcao_m == 0.6
        pasta = "Plots parte 1";
    else
        pasta = "Plots parte 2";
    end

    % Se a pasta ainda não existir, cria
    if ~isfolder(pasta)
        mkdir(pasta)
    end

    % Salva os gráficos como imagens
    saveas(figx, pasta + "\x " + titulo + '.png')
    saveas(figtheta, pasta + "\theta " + titulo + '.png')

end
