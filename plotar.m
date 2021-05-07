function plotar(metodo, t, h, MY, MF, proporcao_m, tracao)

% Plot eixo x
figure;
yyaxis left
plot(t, MY(2, :), '-')
ylabel('Velocidade (m/s)')
yyaxis right
plot(t(2:end), MF(2, :), '-')
ylabel('Aceleração (m/s^2)')
title(tracao + " m " + proporcao_m + " m1 " + (1 - proporcao_m )+ ...
    " " + metodo + " (h = " + h + ")")
grid()
%saveas(gcf,'Barchart.png')

% Plot angulo theta
figure;
yyaxis left
plot(t, MY(3, :), '-')
ylabel('Posição Angular (rad)')
yyaxis right
plot(t, MY(4, :), '-')
ylabel('Velocidade Angular (rad/s)')
title(tracao + " m " + proporcao_m + " m1 " + (1 - proporcao_m )+ ...
    " " + metodo + " (h = " + h + ")")
grid()
%saveas(gcf,'Barchart.png')

end
