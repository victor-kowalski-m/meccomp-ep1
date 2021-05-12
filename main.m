%% Comandos iniciais

clc
clear
close


%% Condições iniciais

Y0 = [0 0 10*pi/180 0];


%% Parte 1

for h=[0.1, 0.01, 0.001]
    
    t = 0:h:20;
    
    for metodo=["euler", "rk2", "rk4"]
        
        sis_eqs = montar_sistema(0.6, "traseira");
        [MY, MF] = resolver(sis_eqs, t, Y0, metodo);
        plotar(metodo, t, MY, MF, 0.6, "traseira")
    
    end
end

%% Parte 2

for tracao=["traseira", "dianteira", "quatro rodas"]
    
    for proporcao_m=[0.2, 0.8]
       
        sis_eqs = montar_sistema(proporcao_m, tracao);
        [MY, MF] = resolver(sis_eqs, t, Y0, "rk4");
        plotar(metodo, t, MY, MF, proporcao_m, tracao)

    end
end