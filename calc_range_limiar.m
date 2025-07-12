%%
% CALCULA UM VETOR DE LIMIARES ADEQUADO PARA CURVAS ROC

%Função: calc_range_limiar

% Autores: André Antônio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberlândia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descrição:
% Calcula um vetor de limiares baseado nas distribuições das
% estatísticas sob as hipóteses H0 e H1, de modo a fornecer thresholds
% apropriados para a construção de curvas ROC ou análise de detecção.

% Entradas:
% N_pts - Número de pontos do vetor de limiares
% Nsimb - Número de símbolos usados na simulação
% Channel - Tipo de canal (ex: 'Rayleigh', 'AWGN', etc.)
% SNRdb - Relação sinal-ruído em dB
% M - Número de subportadoras ou parâmetro específico
% Sepec_Sens_Technique - Técnica de sensoriamento espectral empregada
% limiar - Valor inicial ou auxiliar para cálculo do Pd/Pfa
% nbins - Número de bins usados na estatística (ex: histograma)
% Sample_factor - Fator de amostragem adicional
% Noise_uncertainty_db - Incerteza no nível de ruído em dB
% Error_fc - Erro na frequência central

% Saídas:
% limiar - Vetor (1 x N_pts) contendo os valores dos limiares
% calculados para varredura

function limiar = calc_range_limiar(N_pts,Nsimb,Channel,SNRdb,M,Sepec_Sens_Technique,limiar,nbins,Sample_factor,Noise_uncertainty_db,Error_fc)

display('Calculando range de limiares para ROC');
[~,~,H0,H1] = func_calc_Pd_Pfa_Phase_Dif(1000,Nsimb,Channel,SNRdb,M,Sepec_Sens_Technique,limiar,nbins,Sample_factor,Noise_uncertainty_db,Error_fc);

% Calcula o range total entre o maior valor de H1 e o menor de H0
range_total= max(H1)-min(H0);

% Gera vetor de limiares igualmente espaçados dentro de uma faixa ajustada,
% para evitar extremos onde as curvas ROC podem não ser confiáveis
limiar = linspace((min(H0)+0.05*range_total), (max(H1)-0.05*range_total),N_pts);