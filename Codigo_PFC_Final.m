%  clear all;
%  close all;
%  clc;

%%
% SCRIPT: Spectrum Sensing Phase Difference - ROC Curve

% Autores: André Antônio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberlândia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% DESCRIÇÃO:
% Este script simula e compara o desempenho de técnicas de sensoriamento
% espectral sob condições de ruído impulsivo e incerteza de ruído,
% gerando curvas ROC simuladas para diferentes canais com desvanecimento.
%%

%%
% CONFIGURAÇÕES DE PARÂMETROS

% Número de eventos de Monte Carlo (quantidade de experimentos)
Nmc = 20000;

% Número de amostras coletadas em cada período de sensoriamento
N = 10000;

% SNR média em dB
SNRdb = -10;

% Fator de amostragem fs/Rs
Sample_factor = 100;

% Número de símbolos transmitidos
Nsimb = N/Sample_factor;

% Número de símbolos M-PSK
M = 2;

% Número de BINS para estimar a distribuição de fase
nbins = 40;

% Tipo de canal: 0-Gaussiano; 1-Rayleigh; 2-Rice; 3-Nakagami-m; 4-Alpha-Mu; 5-Kappa-Mu; 6-Eta-Mu
Channel = 1;

% Técnica de sensoriamento espectral: 0-ED, 1-SRD, 2-CRD, 3-PDP-Var, 4-GED
Sepec_Sens_Technique = 0;

% Incerteza de ruído em dB
Noise_uncertainty_db = 0;

% Erro na frequência central
Error_fc = 0;

% Flag para encontrar a faixa válida do limiar
Find_range_limiar = 1;

% Número de pontos para varrer o limiar e construir a curva ROC
N_pts = 500;
%%

%%
% PARÂMETROS DEFAULT DOS CANAIS

CannelParam.Rice.kappa = 3;
CannelParam.Naka.mu = 5;
CannelParam.AlphaMu.alpha = 2;
CannelParam.AlphaMu.mu = 3;
CannelParam.KappaMu.kappa = 10;
CannelParam.KappaMu.mu = 2;
CannelParam.EtaMu.eta = 3;
CannelParam.EtaMu.mu = 0.5;
%%

%%
% PARÂMETROS DO RUÍDO IMPULSIVO (IN)

% Número de RCs
INParam.m = 1;

% Número de amostras coletadas
INParam.n = N;

% Variância do ruído impulsivo = k vezes a variância do ruído térmico
INParam.K = 0;

% Probabilidade de ocorrência do ruído impulsivo
INParam.Pin = 1;

% Probabilidade de cada RC ser afetado pelo ruído impulsivo
INParam.Pcr = 1;

% Amplitude média dos pulsos
INParam.A = 500;

% Variância de cada pulso
INParam.B = 7.5;

% Intervalo médio entre pulsos
INParam.Beta = 100;
%%

%%
%PARÂMETROS A SEREM VARIADOS DOS CANAIS

% Controla qual parâmetro do canal sera variado
% Usa-se 'default' para manter os parâmetros padrão
% Varia canal = 1, parametro = 0
ParametroName = 'default';
parametro = 100000; %valor ignorado no modo default
%%

%%
% EXECUTA A FUNÇÃO PARA GERAR Pd, Pfa E HISTOGRAMAS

[~,~,H_0,H_1] = func_calc_Pd_Pfa_Phase_Dif_vf(Nmc,Nsimb,Channel,SNRdb,M,Sepec_Sens_Technique,0,nbins,Sample_factor,Noise_uncertainty_db,Error_fc,parametro,ParametroName,CannelParam,INParam);
%%

%%
% CONSTRUÇÃO DA CURVA ROC

T_H0 = H_0; % Estatísticas sob hipótese H0 (ausência de sinal)
T_H1 = H_1; % Estatísticas sob hipótese H1 (presença de sinal)

min(T_H0);
max(T_H1);

range_total = max(T_H1)-min(T_H0);

% Definindo faixa de varredura dos limiares com margem de 5% para melhor captura da ROC
% limiar = linspace((min(T_H0)+0.05*range_total), (max(T_H1)-0.05*range_total),Npontos);
limiar = linspace((min(T_H0))-0.05*range_total, (max(T_H1)+0.05*range_total),N_pts);
% limiar = calc_limiar(N_pts,min(T_H0), max(T_H1));

% Calcula Pfa e Pd simulados varrendo os limiares
Pfa = zeros(1,length(limiar));
Pd = zeros(1,length(limiar));

for w = 1: length(limiar)
    Pfa(w) = sum(T_H0>limiar(w))/length(T_H0);
    Pd(w) = sum(T_H1>limiar(w))/length(T_H1);
end
%%

%%
% CÁLCULO DE Pd E Pfa TEÓRICOS

snrlin = 10.^(SNRdb/10);
delta = 2*pi/nbins;
Pd_Pfa_limiar = [Pd', Pfa',limiar'];
Thre = Pd_Pfa_limiar(:,3);

% Pfa teórico para ruído gaussiano (AWGN)
Pfa_teo = qfunc(Thre/sqrt((2*pi-delta)/(4*pi*N)));

if Channel == 0
    % Canal Gaussiano puro (AWGN)
    Pd_teo = qfunc((4*Thre-pi*snrlin)/(4*sqrt(2*pi-delta))*sqrt(4*pi*N));
else
    % Para canais com fading, integra sobre a PDF do canal
    [~,r,fr] = Random_Samples_Fading_model_vf(Channel-1,1,parametro,ParametroName,CannelParam);
    dr = (r(2)-r(1));
    Unit_energy = sum(fr)*dr
    Pd_teo = zeros(1,length(Thre));
    for w = 1: length(Thre)
        Pd_teo(w) = sum(fr.*qfunc((4*Thre(w)-pi*snrlin*r.^2)./(4*sqrt(2*pi-delta))*sqrt(4*pi*N)))*dr;
    end
end

% Cálculo da AUC (Área sob a curva ROC simulada)
AUC_vec = sum(Pd(2:end).*(Pfa(2:end)-Pfa(1:end-1)))*-1
% SNRonSimulation = 10*log10(mean(T_H1)/mean(T_H0)-1)
%%

%%
% HISTOGRAMAS DAS ESTATÍSTICAS SOB H0 E H1 (VISUALIZAÇÃO)

[amp_H0 x_axes_H0] = hist(T_H0,100);
[amp_H1 x_axes_H1] = hist(T_H1,100);
%%

%%
% PLOT DA CURVA ROC (SIMULADA E TEÓRICA)

figure(1);
hold on;
plot(Pfa,Pd,'DisplayName','Sim','Marker','square','LineWidth',2,...
    'LineStyle','--',...
    'Color',[1 0 1]);
xlabel('P_{FA}');
ylabel('P_{D}');
hold on;
grid on;
axis([0 1 0 1]);
plot(Pfa_teo,Pd_teo,'DisplayName','Theo','LineWidth',2,...
    'Color',[0 0 1]);
legend('show')
%%

%%
% GERAÇÃO DOS ARQUIVOS .DAT COM RESULTADOS

if Channel == 0
    Vec_final = [Pd_Pfa_limiar Pfa_teo Pd_teo];
else
    Vec_final = [Pd_Pfa_limiar Pfa_teo Pd_teo'];
end

% scenario = strcat('N',num2str(N),'_SNR_db_',num2str(SNRdb),'_MCe_',num2str(Nmc),'_Channel',num2str(Channel(i)));
% dlmwrite(strcat(scenario,'_Pd_Pfa_Limiar.dat'),Pd_Pfa_Limiar,'delimiter',' ');
% dlmwrite(strcat(scenario,'_Histograma.dat'),Histograma,'delimiter',' ');
%%

%%
% SALVA AUC, GRÁFICO E RESULTADOS PARA ANÁLISE FUTURA

scenario = strcat('N',num2str(N),'_SNR_db_',num2str(SNRdb),'_MCe_',num2str(Nmc));
%legend('show');
dlmwrite(strcat(scenario,'_AUC.dat'),AUC_vec','delimiter',' ');
savefig(strcat(scenario,'Grafico'));
scenario = strcat('N',num2str(N),'_SNR_db_',num2str(SNRdb),'_MCe_',num2str(Nmc),'_Vary_Channel_','All_vecs.dat');
dlmwrite(scenario,Vec_final,'delimiter',' ');
%%

%%
% SALVA CONFIGURAÇÕES DO CENÁRIO PARA DOCUMENTAÇÃO

Param.Channel = Channel;
Param.ParametroName = ParametroName;
Param.Parameter_varying = parametro;
Param.CannelParam = CannelParam;
scenario = strcat('N',num2str(N),'_SNR_db_',num2str(SNRdb),'_MCe_',num2str(Nmc),'_Param.mat')
save(scenario,'Param');
%%