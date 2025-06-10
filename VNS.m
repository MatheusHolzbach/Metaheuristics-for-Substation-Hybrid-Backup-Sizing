% Algoritmo VNS para o dimensionamento de sistemas de backup híbrido
% Universidade Estadual Paulista "Júlio de Mesquita Filho" (UNESP)
% Programa de Pós Graduação em Engenharia Elétrica
% Orientador: Prof. Dr. John Fredy Franco Baquero
% Discente: Matheus Holzbach
% Processo FAPESP n.: 2022/04826-0

clc;
clear;

pasta ='G:\Meu Drive\Mestrado\22_3_Dados\Brasil\';

temperatura = readmatrix("G:\Meu Drive\Mestrado\22_3_Dados\Brasil\temperatura.txt");

wind_speed = readmatrix("G:\Meu Drive\Mestrado\22_3_Dados\Brasil\WindSpeed.txt");

% pasta = 'G:\Meu Drive\Mestrado\22_3_Dados\Canada\';

% temperatura = readmatrix("G:\Meu Drive\Mestrado\22_3_Dados\Canada\temperatura.txt");

% wind_speed = readmatrix("G:\Meu Drive\Mestrado\22_3_Dados\Canada\WindSpeed.txt");

% matriz com os perfis de irradiação (kW/m^2) para 12 meses e 100 perfis com 1440 dados (minuto a minuto)
solar_irrad = zeros(12,100,1440);
for mes = 1 : 12
    for num_profile = 1 : 100
        profile_file = char(strcat(pasta,'PV-',int2str(mes),'-',int2str(num_profile),'.txt'));
        solar_irrad(mes,num_profile,:) = importdata(profile_file);
    end
end

Pcarga = 12;
Tf_media = 5;  % tempo médio de duração da falta
Tf_desvio = 3; % desvio padrão do tempo de duração da falta

N_test = 1;         % número de testes realizados para cada algoritmo
N_Repeticoes = 100; % número de repetições de cada algoritmo por teste

for y = 1 : N_test
    T_max = zeros(1,1);
    Tent_max = zeros(1,1);
    Bs_max = zeros(1,3);
    Bc_max = zeros(1,1);
    Bconf_max = zeros(1,1);
    E_viz_max = zeros(1000,100);
    for j = 1 : N_Repeticoes
        rng('shuffle'); % inicia a semente de gerador de aleatórios
        Solucao_i = randi([30,80],1,3); % gera uma proposta de solução inicial aleatória no intervalo
        Confiabilidade_i = Monte_Carlo(Solucao_i(1), Solucao_i(2), Solucao_i(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed); % avalia a confiabilidade
        Custo_i = round(Valor(Solucao_i(1), Solucao_i(2), Solucao_i(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_i)); % calcula a FO da solução inicial
     
        c = 0; % inicia o contador de iterações(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_a));
        v = 1; % inicia a estrutura de vizinhança
        m = 0; % inicia o contador de melhorias de solução
        w = 0; % inicia o contador de iterações sem melhoria
        RVNS = 10;
        VND = 20;
        tic; % inicia contador de tempo
        while c < 1000 && w < 100 % criterio de parada pelo número máximo de iterações
            if v == 1 % Estratégia de pequena diminuição no número de turbinas eólicas
                Solucao_c = Solucao_i; % restringe o espaço de busca em torno da incumbente
                rng('shuffle');
                M_a = randi([0,RVNS]); % determina um multiplicador aleatório para modificar a solução atual
                Solucao_c(1) = Solucao_c(1) - M_a; % atualiza o valor da nova solução vizinha corrente
                if Solucao_c(1) < 0
                    Solucao_c(1) = 0;
                end
                Confiabilidade_c = Monte_Carlo(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                Custo_c = round(Valor(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_c));
    
                Solucao_a = Solucao_c; % atribui a busca local em torno da solução corente
                Solucao_a(1) = Solucao_a(1) - RVNS;
                if Solucao_c(1) < 0
                    Solucao_c(1) = 0;
                end
                for k = 1 : VND
                    Solucao_a(1) = Solucao_a(1) + 1;
                    Confiabilidade_a = Monte_Carlo(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                    Custo_a = round(Valor(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_a));
                    if Custo_a < Custo_c
                        Solucao_c = Solucao_a; 
                        Custo_c = Custo_a; 
                        Confiabilidade_c = Confiabilidade_a;
                    end
                end
    
                if Custo_c < Custo_i % caso a FO de custo seja melhor que a incumbente e a confiabilidade seja aceitável
                    Solucao_i = Solucao_c; % substituição da solução incumbente em caso de um melhor vizinho encontrado
                    Custo_i = Custo_c; % substituição do custo da solução incumbente
                    Confiabilidade_i = Confiabilidade_c; % substituição da confiabilidade da solução incumbente
                    v = 1; % retorna a primeira estrutura de vizinhança
                    m = m + 1;
                    w = 0;
                else
                    v = 2; % segue para a seguinte estrutura de vizinhança
                    w = w + 1;
                end        
            elseif v == 2 % Estratégia de pequena diminuição no número de baterias
                Solucao_c = Solucao_i; % restringe o espaço de busca em torno da incumbente
                M_a = randi([0,RVNS]); % determina um multiplicador aleatório para modificar a solução atual
                Solucao_c(2) = Solucao_c(2) - M_a; % atualiza o valor da nova solução vizinha corrente
                if Solucao_c(2) < 0
                    Solucao_c(2) = 0;
                end
                Confiabilidade_c = Monte_Carlo(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                Custo_c = round(Valor(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_c));
    
                Solucao_a = Solucao_c; % atribui a busca local em torno da solução corente
                Solucao_a(2) = Solucao_a(2) - RVNS;
                if Solucao_c(2) < 0
                    Solucao_c(2) = 0;
                end
                for k = 1 : VND
                    Solucao_a(2) = Solucao_a(2) + 1;
                    Confiabilidade_a = Monte_Carlo(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                    Custo_a = round(Valor(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_a));
                    if Custo_a < Custo_c
                        Solucao_c = Solucao_a; 
                        Custo_c = Custo_a; 
                        Confiabilidade_c = Confiabilidade_a;
                    end
                end
    
                if Custo_c < Custo_i % caso a FO de custo seja melhor que a incumbente e a confiabilidade seja aceitável
                    Solucao_i = Solucao_c; % substituição da solução incumbente em caso de um melhor vizinho encontrado
                    Custo_i = Custo_c; % substituição do custo da solução incumbente
                    Confiabilidade_i = Confiabilidade_c; % substituição da confiabilidade da solução incumbente
                    v = 1; % retorna para a primeira estrutura de vizinhança
                    m = m + 1;
                    w = 0;
                else
                    v = 3; % segue para a seguinte estrutura de vizinhança
                    w = w + 1;
                end
            elseif v == 3 % Estratégia de pequena diminuição no número de módulos fotovoltaicos
                Solucao_c = Solucao_i; % restringe o espaço de busca em torno da incumbente
                rng('shuffle');
                M_a = randi([0,RVNS]); % determina um multiplicador aleatório para modificar a solução atual
                Solucao_c(3) = Solucao_c(3) - M_a; % atualiza o valor da nova solução vizinha corrente
                if Solucao_c(3) < 0
                    Solucao_c(3) = 0;
                end
                Confiabilidade_c = Monte_Carlo(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                Custo_c = round(Valor(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_c));
    
                Solucao_a = Solucao_c; % atribui a busca local em torno da solução corente
                Solucao_a(3) = Solucao_a(3) - RVNS;
                if Solucao_c(3) < 0
                    Solucao_c(3) = 0;
                end
                for k = 1 : VND
                    Solucao_a(3) = Solucao_a(3) + 1;
                    Confiabilidade_a = Monte_Carlo(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                    Custo_a = round(Valor(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_a));
                    if Custo_a < Custo_c
                        Solucao_c = Solucao_a; 
                        Custo_c = Custo_a; 
                        Confiabilidade_c = Confiabilidade_a;
                    end
                end
    
                if Custo_c < Custo_i % caso a FO de custo seja melhor que a incumbente e a confiabilidade seja aceitável
                    Solucao_i = Solucao_c; % substituição da solução incumbente em caso de um melhor vizinho encontrado
                    Custo_i = Custo_c; % substituição do custo da solução incumbente
                    Confiabilidade_i = Confiabilidade_c; % substituição da confiabilidade da solução incumbente
                    v = 1; % retorna para a primeira estrutura de vizinhança
                    m = m + 1;
                    w = 0;
                else
                    v = 4; % segue para a seguinte estrutura de vizinhança
                    w = w + 1;
                end
            elseif v == 4 % Estratégia de pequeno aumento no número de módulos fotovoltaicos
                Solucao_c = Solucao_i; % restringe o espaço de busca em torno da incumbente
                rng('shuffle');
                M_a = randi([0,RVNS]); % determina um multiplicador aleatório para modificar a solução atual
                Solucao_c(3) = Solucao_c(3) + M_a; % atualiza o valor da nova solução vizinha corrente
                
                Confiabilidade_c = Monte_Carlo(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                Custo_c = round(Valor(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_c));
    
                Solucao_a = Solucao_c; % atribui a busca local em torno da solução corente
                Solucao_a(3) = Solucao_a(3) - RVNS;
                if Solucao_c(3) < 0
                    Solucao_c(3) = 0;
                end
                for k = 1 : VND
                    Solucao_a(3) = Solucao_a(3) + 1;
                    Confiabilidade_a = Monte_Carlo(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                    Custo_a = round(Valor(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_a));
                    if Custo_a < Custo_c
                        Solucao_c = Solucao_a; 
                        Custo_c = Custo_a; 
                        Confiabilidade_c = Confiabilidade_a;
                    end
                end
    
                if Custo_c < Custo_i % caso a FO de custo seja melhor que a incumbente e a confiabilidade seja aceitável
                    Solucao_i = Solucao_c; % substituição da solução incumbente em caso de um melhor vizinho encontrado
                    Custo_i = Custo_c; % substituição do custo da solução incumbente
                    Confiabilidade_i = Confiabilidade_c; % substituição da confiabilidade da solução incumbente
                    v = 1; % retorna para a primeira estrutura de vizinhança
                    m = m + 1;
                    w = 0;
                else
                    v = 5; % retora para a primeira estrutura de vizinhança
                    w = w + 1;
                end
            elseif v == 5 % Estratégia de pequeno incremento no número de geradores eólicos
                Solucao_c = Solucao_i; % restringe o espaço de busca em torno da incumbente
                rng('shuffle');
                M_a = randi([0,RVNS]); % determina um multiplicador aleatório para modificar a solução atual
                Solucao_c(1) = Solucao_c(1) + M_a; % atualiza o valor da nova solução vizinha corrente
    
                Confiabilidade_c = Monte_Carlo(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                Custo_c = round(Valor(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_c));
    
                Solucao_a = Solucao_c; % atribui a busca local em torno da solução corente
                Solucao_a(1) = Solucao_a(1) - RVNS;
                if Solucao_c(1) < 0
                    Solucao_c(1) = 0;
                end
                for k = 1 : VND
                    Solucao_a(1) = Solucao_a(1) + 1;
                    Confiabilidade_a = Monte_Carlo(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                    Custo_a = round(Valor(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_a));
                    if Custo_a < Custo_c
                        Solucao_c = Solucao_a; 
                        Custo_c = Custo_a; 
                        Confiabilidade_c = Confiabilidade_a;
                    end
                end
    
                if Custo_c < Custo_i % caso a FO de custo seja melhor que a incumbente e a confiabilidade seja aceitável
                    Solucao_i = Solucao_c; % substituição da solução incumbente em caso de um melhor vizinho encontrado
                    Custo_i = Custo_c; % substituição do custo da solução incumbente
                    Confiabilidade_i = Confiabilidade_c; % substituição da confiabilidade da solução incumbente
                    v = 1; % retorna para a primeira estrutura de vizinhança
                    m = m + 1;
                    w = 0;
                else
                    v = 6; % segue para a seguinte estrutura de vizinhança
                    w = w + 1;
                end
            elseif v == 6 % Estratégia de pequeno incremento no número de baterias
                Solucao_c = Solucao_i; % restringe o espaço de busca em torno da incumbente
                rng('shuffle');
                M_a = randi([0,RVNS]); % determina um multiplicador aleatório para modificar a solução atual
                Solucao_c(2) = Solucao_c(2) + M_a; % atualiza o valor da nova solução vizinha corrente
    
                Confiabilidade_c = Monte_Carlo(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                Custo_c = round(Valor(Solucao_c(1), Solucao_c(2), Solucao_c(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_c));
    
                Solucao_a = Solucao_c; % atribui a busca local em torno da solução corente
                Solucao_a(2) = Solucao_a(2) - RVNS;
                if Solucao_c(2) < 0
                    Solucao_c(2) = 0;
                end
                for k = 1 : VND
                    Solucao_a(2) = Solucao_a(2) + 1;
                    Confiabilidade_a = Monte_Carlo(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed);
                    Custo_a = round(Valor(Solucao_a(1), Solucao_a(2), Solucao_a(3), Pcarga, solar_irrad, wind_speed, Confiabilidade_a));
                    if Custo_a < Custo_c
                        Solucao_c = Solucao_a; 
                        Custo_c = Custo_a; 
                        Confiabilidade_c = Confiabilidade_a;
                    end
                end
    
                if Custo_c < Custo_i % caso a FO de custo seja melhor que a incumbente e a confiabilidade seja aceitável
                    Solucao_i = Solucao_c; % substituição da solução incumbente em caso de um melhor vizinho encontrado
                    Custo_i = Custo_c; % substituição do custo da solução incumbente
                    Confiabilidade_i = Confiabilidade_c; % substituição da confiabilidade da solução incumbente
                    v = 1; % retorna para a primeira estrutura de vizinhança
                    m = m + 1;
                    w = 0;
                else
                    v = 1; % segue para a seguinte estrutura de vizinhança
                    w = w + 1;
                end
            end
            c = c + 1;
        end
        toc; % finaliza contador de tempo
        T_max(j) = toc;
        Tent_max(j) = c;
        Bs_max(j,:) = Solucao_i;
        Bc_max(j) = Custo_i;
        Bconf_max(j) = Confiabilidade_i;
        fprintf("Verificação: Iteração %0.0f/%0.0f VNS %0.0f \n", j, N_Repeticoes, y);
    end
    
    n_arq1 = char(strcat('T_max_VNS_',int2str(y),'.txt'));
    n_arq2 = char(strcat('Tent_max_VNS_',int2str(y),'.txt'));
    n_arq3 = char(strcat('Bs_max_VNS_',int2str(y),'.txt'));
    n_arq4 = char(strcat('Bc_max_VNS_',int2str(y),'.txt'));
    n_arq5 = char(strcat('Bconf_max_VNS_',int2str(y),'.txt'));

    writematrix(T_max,n_arq1);
    writematrix(Tent_max,n_arq2);
    writematrix(Bs_max,n_arq3);
    writematrix(Bc_max,n_arq4);
    writematrix(Bconf_max,n_arq5);
end

function [Confiabilidade] = Monte_Carlo(NWT, NB, NPV, Pcarga, Tf_media, Tf_desvio, solar_irrad, temperatura, wind_speed)
    rng(4); % inicia a semente de gerador de aleatórios
    Noct = 44;
    Voc = 65.8;
    Isc = 4.89;
    Kv = -0.174;
    Ki = 0.00182;
    Vmppt = 58.0;
    Imppt = 5.70;
    FF = Vmppt * Imppt / (Voc * Isc);

    VCut_in = 2; % Velocidade do vento de corte inferior
    VCut_on = 20; % Velocidade do vento de corte superior
    Vr = 9; % Velocidade nominal do vento
    PWT_max = 2; % Potência máxima da turbina eólica
    DT = 12; % Diâmetro da turbina
    D_ar = 1.2754; % Densidade relativa média do ar

    CP_B = 2; % capacidade da Bateria (kW)
    DoD = 0.9; % Deep-of-Discharge da bateria
    Ef_Desc_B = 0.95; % eficiência de descarga da bateria

    T_Indisp_total = 0;
    Tf_d_total = 0;
    Nf_s = 500;
    for Nf = 1 : Nf_s
        for Mf = 1 : 12
            Tf_d = round(max(0,Tf_media + Tf_desvio*randn(1,1)) * 60); % tempo de duração da falta aleatório
            if Tf_d >= 1440
                Tf_d = 1439;
            end
            Tf_i = randi(1440 - Tf_d); % tempo de inicio da falta aleatório
            Tf_f = Tf_i + Tf_d; % tempo da finalização da falta
            Cf = randi(100);
            Df = Mf*30 - floor(rand()*30);
    
            E_B = Ef_Desc_B * DoD * NB * CP_B; % energia disponível na bateria
            E_Disp = E_B; % energia inicial disponível no sistema híbrido
    
            Tf_s = Tf_i;
            Autonomia = 1;
            while Tf_s < Tf_f && Autonomia == 1
                Irrad_solar = solar_irrad(Mf,Cf,Tf_s);
                Ta = temperatura((Df-1) * 24 + ceil((Tf_i+Tf_s)/60));
                Tc = Ta + (Noct - 20) / 0.8 * Irrad_solar;
                Ic = Irrad_solar * (Isc + Ki * (Tc - 25) );
                Vc = Voc + Kv * Tc;
                E_PV = ((NPV * FF * Vc * Ic) / 60) / 1000;
    
                VW = wind_speed((Df-1) * 24 + ceil((Tf_i+Tf_s)/60));
                PWT_r = PWT_max * 0.5 * D_ar * ((pi * DT.^2) ./ 4) * Vr.^3;
                PWT = VW.^3 * ((PWT_r) ./ (Vr.^3 - VCut_in.^3)) - PWT_r * ((VCut_in.^3) ./ (Vr.^3 - VCut_in.^3));
                if VW < VCut_in
                    E_WT = 0;
                elseif VW >= VCut_in && VW < Vr
                    E_WT = NWT * (PWT ./1000 ./60);
                elseif VW >= Vr && VW <= VCut_on
                    E_WT = NWT * (PWT_r ./1000 ./60);
                elseif VW > VCut_on
                    E_WT = 0;
                end
    
                E_Disp = E_Disp + E_PV + E_WT - Pcarga/60;
                Tf_s = Tf_s + 1;
                if E_Disp <= 0
                    Autonomia = 0;
                end
            end
            if E_Disp <= 0
                Tf_s = Tf_s - 1; % Corrige o incremento de tempo de simulação
            end
            T_Indisp = (Tf_f - Tf_s)/60; % tempo de indisponibilidade em horas
            T_Indisp_total = T_Indisp_total + T_Indisp; % atualiza o tempo total indisponível
            Tf_d_total = Tf_d_total + (Tf_d / 60);
        end
    end
    Indice = (T_Indisp_total / Tf_d_total) * 100; % calcula o índice (%) como a razão do tempo total indisponível para o número de horas totais dos anos de simulação
    Confiabilidade = 100 - Indice;
end

function [VLP] = Valor(NWT, NB, NPV, Pcarga, solar_irrad, wind_speed, Conf)
    cus_bat = 420; % USD/kWh
    cus_pv = 312; % USD/unidade
    cus_wt = 1500; % USD/kW
    pot_pv = 0.33;
    pot_wt = 1.5;
    cus_inv_pv = 105;
    cus_inv_wt = 140;
    man_bat = 0.015;
    man_pv = 0.01;
    man_wt = 0.015;
    man_inv_pv = 0.015;
    man_inv_wt = 0.015;
    vida_util = 20; % anos
    mod_bat = 2; % tamanho do módulo da bateria em kWh
    C_DG = 0.03201;
    juros = 0.06;
    energy_price = 0.05;    % preço da energia vendida (alguns centavos de dólar)
    fator_vp = 0;
    for k = 1 : vida_util
        fator_vp = fator_vp + 1/(1+juros)^k;
    end

    lucroPV = sum(sum(sum(solar_irrad)))/12/100/60*365*pot_pv*energy_price*fator_vp;
    lucroWT = sum(wind_speed)/24/365*pot_wt*energy_price*fator_vp;
    O_DG = (100-Conf)*8760*Pcarga*C_DG*fator_vp;

    p_inv_pv = max (NPV*pot_pv,Pcarga);
    p_inv_wt = max (NWT*pot_wt,Pcarga);
    custo = NB*mod_bat*cus_bat + NPV*cus_pv + p_inv_pv*cus_inv_pv + NWT*pot_wt*cus_wt + p_inv_wt*cus_inv_wt;
        
    manutencao = 0;
    for k = 1 : vida_util
        manutencao = manutencao + (NB*mod_bat*man_bat + NPV*man_pv + man_inv_pv*p_inv_pv + NWT*pot_wt*man_wt + man_inv_wt*p_inv_wt)/(1+juros)^k;
    end
    VLP = (custo + manutencao + O_DG - lucroPV*NPV - lucroWT*NWT);
end