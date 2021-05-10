%% Comandos iniciais

clc
clear
close

%% Globais

% Constantes
mtotal = 1939;
L = 2.95;
I = 1;
w = 10;
mi = 0.42;
beta = 0.02;
g = 9.8;

% Condições iniciais
Y0 = [0 0 10*pi/180 0];

save('globais.mat')

%% Parte 1

for h=[0.1, 0.01, 0.001]
    
    t = 0:h:20;
    
    for metodo=["euler", "rk2", "rk4"]
        
        resolver(metodo, h, t, "traseira", 0.6)
    
    end
end

%% Parte 2

for tracao=["traseira", "dianteira", "quatro rodas"]
    
    for proporcao_m=[0.2, 0.8]

        resolver("rk4", h, t, tracao, proporcao_m)

    end
end