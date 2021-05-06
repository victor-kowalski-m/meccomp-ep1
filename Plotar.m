function Plotar(metodo, t, h, MY, MF)

figure;
yyaxis left
plot(t, MY(2, :), '-')
ylabel('Velocidade (m/s)')
yyaxis right
plot(t(2:end), MF(2, :), '-')
ylabel('Aceleração (m/s^2)')
title(metodo + " Posição e Velocidade (h = " + h + ")")
grid()

figure;
yyaxis left
plot(t, MY(3, :), '-')
ylabel('Posição Angular (rad)')
yyaxis right
plot(t, MY(4, :), '-')
ylabel('Velocidade Angular (rad/s)')
title(metodo + " Posição e Velocidade Angulares (h = " + h + ")")
grid()

end