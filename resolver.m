function resolver(metodo, h, t, tracao, proporcao_m)

% Carrega variáveis globais
load('globais.mat')
load('globais.mat', 'beta')

% Distribuição de massas
m = proporcao_m*mtotal;
m1 = (1 - proporcao_m)*mtotal;

% Forças atuantes nos pneus
F = beta *m*g;
F1 = mi*m1*g;

% Matrizes iniciais
MY = [];
MY(:, 1) = Y0;
MF = [];

% Define numero de iteracoes
iteracoes = length(t)-1

% Resolve conforme o método escohido
if metodo == "euler"
    for i=1:iteracoes

        Y = MY(:,i);

        K = [
            Y(2);
            -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
            + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
            Y(4);
            -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
            - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
            /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
            ];

        Y = MY(:,i) + h * K;

        MY(:,i+1) = Y;
        MF(:,i) = K;
        
    end

elseif metodo == "rk2"
    for i=1:iteracoes

        Y = MY(:,i);

        K1 = [
            Y(2);
            -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
            + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
            Y(4);
            -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
            - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
            /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
            ];

        Y = MY(:,i) + h/2 * K1;

        K2 = [
            Y(2);
            -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
            + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
            Y(4);
            -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
            - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
            /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
            ];

        Y = MY(:,i) + h * K2;
        MY(:,i+1) = Y;
        MF(:,i) = K2;
        
    end

elseif metodo == "rk4"
    for i = 1:iteracoes

        Y = MY(:,i);

        K1 = [
            Y(2);
            -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
            + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
            Y(4);
            -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
            - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
            /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
            ];

        Y = MY(:,i) + h/2 * K1;

        K2 = [
            Y(2);
            -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
            + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
            Y(4);
            -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
            - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
            /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
            ];

        Y = MY(:,i) + h/2 * K2;

        K3 = [
            Y(2);
            -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
            + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
            Y(4);
            -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
            - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
            /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
            ];

        Y = MY(:,i) + h * K3;

        K4 = [
            Y(2);
            -(F*L - F1*L*cos(Y(3)) - 2*I*w*sin(Y(3))*Y(4) ...
            + L^2*m1*cos(Y(3))*Y(4)^2)/(L*(mtotal - m1*sin(Y(3))^2));
            Y(4);
            -((F*L*m1*sin(Y(3)) - F1*L*m1*cos(Y(3))*sin(Y(3)) ...
            - 2*I*mtotal*w*Y(4) + L^2*m1^2*cos(Y(3))*sin(Y(3))*Y(4)^2) ...
            /(L^2*m1*(-mtotal + m1*sin(Y(3))^2)))
            ];

        Y = MY(:,i) + h/6 * (K1 + 2*K2 + 2*K3 + K4);
        MY(:, i+ 1) = Y;
        MF(:, i) = 1/6 * (K1 + 2*K2 + 2*K3 + K4);

    end
    
else
    disp("Metodo inválido.")    
end

% Plota e salva resultados
plotar(metodo, t, h, MY, MF, proporcao_m, tracao)

end