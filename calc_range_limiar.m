%%
% CALCULA UM VETOR DE LIMIARES ADEQUADO PARA CURVAS ROC

%Fun��o: calc_range_limiar

% Autores: Andr� Ant�nio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberl�ndia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descri��o:
% Calcula um vetor de limiares baseado nas distribui��es das
% estat�sticas sob as hip�teses H0 e H1, de modo a fornecer thresholds
% apropriados para a constru��o de curvas ROC ou an�lise de detec��o.

% Entradas:
% N_pts - N�mero de pontos do vetor de limiares
% Nsimb - N�mero de s�mbolos usados na simula��o
% Channel - Tipo de canal (ex: 'Rayleigh', 'AWGN', etc.)
% SNRdb - Rela��o sinal-ru�do em dB
% M - N�mero de subportadoras ou par�metro espec�fico
% Sepec_Sens_Technique - T�cnica de sensoriamento espectral empregada
% limiar - Valor inicial ou auxiliar para c�lculo do Pd/Pfa
% nbins - N�mero de bins usados na estat�stica (ex: histograma)
% Sample_factor - Fator de amostragem adicional
% Noise_uncertainty_db - Incerteza no n�vel de ru�do em dB
% Error_fc - Erro na frequ�ncia central

% Sa�das:
% limiar - Vetor (1 x N_pts) contendo os valores dos limiares
% calculados para varredura

function limiar = calc_range_limiar(N_pts,Nsimb,Channel,SNRdb,M,Sepec_Sens_Technique,limiar,nbins,Sample_factor,Noise_uncertainty_db,Error_fc)

display('Calculando range de limiares para ROC');
[~,~,H0,H1] = func_calc_Pd_Pfa_Phase_Dif(1000,Nsimb,Channel,SNRdb,M,Sepec_Sens_Technique,limiar,nbins,Sample_factor,Noise_uncertainty_db,Error_fc);

% Calcula o range total entre o maior valor de H1 e o menor de H0
range_total= max(H1)-min(H0);

% Gera vetor de limiares igualmente espa�ados dentro de uma faixa ajustada,
% para evitar extremos onde as curvas ROC podem n�o ser confi�veis
limiar = linspace((min(H0)+0.05*range_total), (max(H1)-0.05*range_total),N_pts);