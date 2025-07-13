%%
% CALCULA Pd E Pfa UTILIZANDO TÉCNICAS DE SENSORIAMENTO POR DIFERENÇA DE FASE

% Função: func_calc_Pd_Pfa_Phase_Dif_vf

% Autores: André Antônio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberlândia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descrição:
% Calcula as probabilidades de detecção (Pd) e falso alarme (Pfa)
% para diferentes técnicas de sensoriamento espectral baseadas em sinais,
% utilizando estatísticas extraídas da diferença de fase entre amostras do
% sinal recebido. O cálculo é realizado por meio de simulações de Monte Carlo,
% considerando diferentes cenários de canal, ruído e técnicas.

% Entradas:
% Nmc - Número de eventos de Monte Carlo
% Nsimb - Número de símbolos transmitidos
% Channel - Tipo de canal (0: Gaussiano, >=1: Rayleigh lento + IN)
% SNRdb - SNR em dB
% M - Ordem da modulação M-PSK
% Sepec_Sens_Technique - Técnica de sensoriamento (0: energia, 1: SRD, 2: CRD, 3: CRD+energia)
% Thre - Limiar de detecção
% nbins - Número de bins para estimar a pdf da diferença de fase
% Sample_factor - Fator de superamostragem
% Noise_uncertainty_db - Incerteza de ruído em dB
% Error_fc - Erro na frequência central
% parametro, ParametroName, CannelParam - Parâmetros do canal
% INParam - Struct com parâmetros do ruído impulsivo

% Saídas:
% Pd - Probabilidade de detecção
% Pfa - Probabilidade de falso alarme
% T_H0 - Vetor das estatísticas sob H0 (canal ocioso)
% T_H1 - Vetor das estatísticas sob H1 (canal ocupado)
%%

function [Pd,Pfa,T_H0,T_H1] = func_calc_Pd_Pfa_Phase_Dif_vf(Nmc,Nsimb,Channel,SNRdb,M,Sepec_Sens_Technique,Thre,nbins,Sample_factor,Noise_uncertainty_db,Error_fc,parametro,ParametroName,CannelParam,INParam)

%%
% INICIALIZA VETORES PARA ARMAZENAR ESTATÍSTICAS SOB H0 (CANAL OCIOSO) E H1 (CANAL OCUPADO)

T_H0 = zeros(1,Nmc/2);
T_H1 = zeros(1,Nmc/2);
%%

%%
% CONFIGURAÇÕES DO SISTEMA

% Eventos de Montecarlo
% Nmc = 1000;

% Número de símbolos transmitidos
% Nsimb = 1000;

% Canal Rayleigh lento (1) ou gaussino (0)
% Channel = 1;

% SNR na recepção
% SNRdb = -20;

% Ordem da modulação MPSK considerada
% M = 2;
%%

%%
% DEFINIÇÕES DAS TAXAS ENVOLVIDAS

% Taxa de bit
Rb = 1;

% Taxa de símbolo
Rs = Rb ./ log2(M);

% Tempo de símbolo
Tsimb = 1/Rs;

% Fator de amostragem (quantas vezes a frequência de amostragem é maior que a taxa de símbolo)
K = Sample_factor;

% Taxa de amostragem
Fs = K*Rs;

% Tempo de amostragem
Ts = 1/Fs;

% Frequência do oscilador local
fosc = Fs/4;
%%

%%
% INICIALIZA CONTADORES PARA PD E PFA

contpd = 0;
contpfa = 0;
T_H0 = zeros(1,Nmc/2);
T_H1 = zeros(1,Nmc/2);
%%

%%
% INCERTEZA DE RUÍDO EM ESCALA LINEAR

Noise_uncertainty = 10.^(Noise_uncertainty_db/10);
%%

%%
% LOOP MONTE CARLO

for i = 1 : Nmc
    
    %%
    % DEFINE SE É HIPÓTESE H0 (CANAL OCIOSO) OU H1 (CANAL OCUPADO)
    
    if i <= Nmc/2
        TxEnable = 0;
    else
        TxEnable = 1;
    end
    %%
    
    %%
    % GERAÇÃO DOS SÍMBOLOS TRANSMITIDOS
    
    Data = Gen_Simb(M,Nsimb);
    % if Channel == 1
    % Data = Data .*(((randn(1,Nsimb)+1i*randn(1,Nsimb))*sig_ray));
    % Data = Data .*raylrnd(1,1,Nsimb);
    % end
    % if Channel == 1
    % h = (1/sqrt(2)*(randn(1,Nsimb)+1i*randn(1,Nsimb)));
    % else
    % h = ones(1,Nsimb);
    % end
    %%
    
    %%
    % GERAÇÃO DOS DADOS AMOSTRADOS
    
    Data_Kx = zeros(1,K*Nsimb);
    %h_Kx = zeros(1,K*Nsimb);
    
    for j = 1 : K
        Data_Kx(j:K:end)=Data;
        %h_Kx(j:K:end)= h;
    end
    
    % Data_Kx = repelem(Data,K);
    % h_Kx = repelem(h,K);
    % if Channel == 1
    % Data_Kx = Data_Kx .*(((randn(1,length(Data_Kx))+1i*randn(1,length(Data_Kx)))*sig_ray));
    % Data = Data .*raylrnd(1,1,Nsimb);
    % end
    %%
    
    %%
    % CRIAÇÃO DO VETOR DE TEMPO PARA O SINAL MODULADO
    
    t = [0: 1/Fs : (length(Data_Kx) - 1)*1/Fs];
    % cosseno = sqrt(2/Tsimb)*cos(2*pi*fosc.*t);
    % seno = sqrt(2/Tsimb)*sin(2*pi*fosc.*t);
    %%
    
    %%
    % MODULAÇÃO COMPLEXA COM O OSCILADOR LOCAL (TRANSLAÇÃO EM FREQUÊNCIA)
    
    % Data_Mod = real(Data_Kx) .*cosseno + 1i*imag(Data_Kx).*seno;
    
    % Símbolo multiplicado pela exponencial complexa e^(j2*pi*fc/fst)
    % Data_Mod = Data_Kx.*(cosseno+1i*seno);
    % Data_Mod = Data_Kx.*exp(1i*2*pi*fosc.*t);
    % Data_Mod = abs(Data_Kx).*exp(1i*(2*pi*fosc.*t+phase(Data_Kx))).*(h_Kx);
    Data_Mod = Data_Kx.*exp(1i*2*pi*fosc.*t);%.*(h_Kx);
    % Data_Mod = Data_Kx.*(h_Kx);
    % if Channel == 1
    % h = [1 0 0 0.5 2];
    % Data_Mod = filter(h,1,Data_Mod);
    % end
    % Pot_sinal = sum(real(Data_Mod).^2 + imag(Data_Mod).^2)./length(Data_Mod);
    % if Channel == 1
    % Pot_sinal = 2*sig_ray.^2;
    % end
    % Pot_sinal = 1;
    % Noise = sqrt(1/((10.^(SNRdb/10))/Pot_sinal))*(1/sqrt(2)*randn(1,length(Data_Mod)) + 1i* 1/sqrt(2)*randn(1,length(Data_Mod)));
    % Data_rx = TxEnable*Data_Mod + Noise;
    %%
    
    %%
    % NORMALIZAÇÃO DE ENERGIA E ADIÇÃO DE RUÍDO AWGN COM INCERTEZA DE RUÍDO
    
    Pot_sinal = sum(real(Data_Mod).^2 + imag(Data_Mod).^2)./length(Data_Mod);
    Data_Mod = Data_Mod/sqrt(Pot_sinal)*sqrt(10.^(SNRdb/10));
    Noise = sqrt(unifrnd(1/Noise_uncertainty,Noise_uncertainty))*(1/sqrt(2)*randn(1,length(Data_Mod)) + 1i* 1/sqrt(2)*randn(1,length(Data_Mod)));
    % Data_rx = TxEnable*Data_Mod + Noise;
    
    % Aplicação do canal (fading) e ruído impulsivo se necessário
    if Channel >= 1
        % Data_rx = TxEnable*Data_Mod*(1/sqrt(2)*(randn(1,1)+1i*randn(1,1))) + Noise;
        [x,~,~] = Random_Samples_Fading_model_vf(Channel-1,1,parametro,ParametroName,CannelParam);
        Impusive_Noise = Gen_Impulsive_Noise(INParam.m,INParam.n,0,0,INParam.K,INParam.Pin,INParam.Pcr,INParam.A,INParam.B,INParam.Beta,2,0);
        Data_rx = TxEnable*Data_Mod*(x) + Noise + Impusive_Noise;
    else
        Impusive_Noise = Gen_Impulsive_Noise(INParam.m,INParam.n,0,0,INParam.K,INParam.Pin,INParam.Pcr,INParam.A,INParam.B,INParam.Beta,2,0);
        Data_rx = TxEnable*Data_Mod + Noise + Impusive_Noise;
    end
    %%
    
    %%
    % ESTIMAÇÃO DAS FASES INSTANTÂNEAS
    
    %fi_n = mod(phase(Data_rx),2*pi); % Estava demorando muito
    %fi_n = mod(atan2(imag(Data_rx),real(Data_rx)),2*pi);
    fi_n = atan2(imag(Data_rx),real(Data_rx));
    %%
    
    %%
    % CALCULA A DIFERENÇA DE FASE ENTRE AMOSTRAS CONSECUTIVAS
    
    theta_n = mod((fi_n(2:end)-fi_n(1:end-1)),2*pi);
    
    %Plota histogramas da fase do sinal recebido e da diferença de fase do sinal recebido
    % figure; hist(fi_n,50);legend('pdf de phi_n');
    % figure; hist(theta_n,50);legend('pdf de theta');
    %%
    
    %%
    % CALCULA A PDF EMPÍRICA DA DIFERENÇA DE FASE
    
    [xpdf_empirical,ypdf_empirical] = fcn_coder_pdf(theta_n, nbins);
    % figure; stem(xpdf_empirical,ypdf_empirical,'--*b')
    % var(ypdf_empirical);
    %%
    
    %%
    % PLOTANDO AS PDFS TEÓRICAS
    
    % theta = 0:1/100*pi:2*pi;
    % snrdb = SNRdb;
    % snr_lin = 10.^(snrdb/10);
    % fc_div_fs = fosc/Fs;
    % f_theta_eq10 = 1/(2*pi)+ (snr_lin)/4 *cos(2*pi*fc_div_fs - theta)+ (snr_lin).^2/(2*pi*(snr_lin+1)) *cos(2*(2*pi*fc_div_fs - theta));
    % f_theta_eq11 = 1/(2*pi)+ (snr_lin)/4 *cos(2*pi*fc_div_fs - theta);
    % hold on
    % plot(theta,f_theta_eq10,'+r');
    % plot(theta,f_theta_eq11,'--r');
    % legend('pdf empirical','pdf theoretical eq10','pdf theoretical eq11')
    % xlabel('Diferença de fase [rad]')
    % ylabel('Pr')
    %%
    
    %%
    % CÁLCULO DA VARIÁVEL DE DECISÃO COM BASE NA TÉCNICA DE DETECÇÃO ESCOLHIDA
    
    [fmax, pxmax]= max(ypdf_empirical);
    
    % Detecção de energia
    if Sepec_Sens_Technique == 0
        T = sum(abs(Data_rx).^2);
        
        % SRD baseado nas tangentes das bordas da PDF suavizada
    elseif Sepec_Sens_Technique == 1
        ypdf_empirical = smooth(xpdf_empirical,ypdf_empirical,0.86,'loess');
        tan_am = (ypdf_empirical(floor(nbins/2)+1) - ypdf_empirical(1))/xpdf_empirical(floor(nbins/2)+1);
        tan_bm = (ypdf_empirical(floor(nbins/2)+1) - ypdf_empirical(nbins))/(xpdf_empirical(floor(nbins/2)+1)- 2*pi);
        % tan_am = (fmax - ypdf_empirical(1))/xpdf_empirical(pxmax);
        % tan_bm = (fmax - ypdf_empirical(nbins))/(xpdf_empirical(floor(nbins/2)+1)- 2*pi);
        T = ((tan_am)-(tan_bm))/2;
        
        % CRD
    elseif Sepec_Sens_Technique == 2
        % thetai = linspace(0,2*pi-2*pi/nbins,nbins);
        thetai = xpdf_empirical;
        T = sum(ypdf_empirical.*cos(2*pi*(1+Error_fc)*fosc/Fs - thetai)*2*pi/nbins);
        
        % Variante CRD com energia do sinal adicionada
    elseif Sepec_Sens_Technique == 3
        % ypdf_empirical = smooth(xpdf_empirical,ypdf_empirical,0.86,'loess')';
        % T = var((ypdf_empirical).^2); %PDP-Variance
        % thetai = xpdf_empirical;
        % T = sum(abs(ypdf_empirical-mean(ypdf_empirical)).^5);
        % T = var((ypdf_empirical).^2) /length(ypdf_empirical);
        % Esse do NMSE deu um bom resultado
        % T = sum((ypdf_empirical-cos(2*pi*(1+Error_fc)*fosc/Fs - thetai)).^2)/sum((cos(2*pi*(1+Error_fc)*fosc/Fs - thetai)).^2);
        % T = 1/T;
        % T = sum(((ypdf_empirical-(2*pi)/nbins)).^2);
        % T = 1/T
        thetai = xpdf_empirical;
        Energy = sum(abs(Data_rx).^2);
        T = sum((ypdf_empirical).*cos(2*pi*(1+Error_fc)*fosc/Fs - thetai)*2*pi/nbins)+Energy;
        
        % Fallback: detecção com potência fracionária
    else
        T = sum(abs(Data_rx).^1.1);
    end
    %%
    
    %%
    % ENCONTRANDO A THRE
    
    % Thre = 0.1127/sqrt(N)*(qfuncinv(Pfa))
    %%
    
    %%
    % ARMAZENA ESTATÍSTICA E CONTA PD OU PFA
    
    if i<= Nmc/2
        T_H0(i) = T;
        
        if T  >= Thre
            
            % Display('Ocupado');
            contpfa = contpfa +1;
        else
            % Display('vazio');
        end
    else
        T_H1(i-Nmc/2) = T;
        if T  >= Thre
            % Display('Ocupado');
            contpd = contpd +1;
        else
            % Display('Vazio');
        end
    end
end
%%

%%
% CALCULA AS PROBABLIDADES PD E PFA

Pfa = contpfa/(Nmc/2);
Pd = contpd/(Nmc/2);
%%
