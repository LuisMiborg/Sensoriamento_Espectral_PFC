%%
% CALCULA Pd E Pfa UTILIZANDO T�CNICAS DE SENSORIAMENTO POR DIFEREN�A DE FASE

% Fun��o: func_calc_Pd_Pfa_Phase_Dif_vf

% Autores: Andr� Ant�nio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberl�ndia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descri��o:
% Calcula as probabilidades de detec��o (Pd) e falso alarme (Pfa)
% para diferentes t�cnicas de sensoriamento espectral baseadas em sinais,
% utilizando estat�sticas extra�das da diferen�a de fase entre amostras do
% sinal recebido. O c�lculo � realizado por meio de simula��es de Monte Carlo,
% considerando diferentes cen�rios de canal, ru�do e t�cnicas.

% Entradas:
% Nmc - N�mero de eventos de Monte Carlo
% Nsimb - N�mero de s�mbolos transmitidos
% Channel - Tipo de canal (0: Gaussiano, >=1: Rayleigh lento + IN)
% SNRdb - SNR em dB
% M - Ordem da modula��o M-PSK
% Sepec_Sens_Technique - T�cnica de sensoriamento (0: energia, 1: SRD, 2: CRD, 3: CRD+energia)
% Thre - Limiar de detec��o
% nbins - N�mero de bins para estimar a pdf da diferen�a de fase
% Sample_factor - Fator de superamostragem
% Noise_uncertainty_db - Incerteza de ru�do em dB
% Error_fc - Erro na frequ�ncia central
% parametro, ParametroName, CannelParam - Par�metros do canal
% INParam - Struct com par�metros do ru�do impulsivo

% Sa�das:
% Pd - Probabilidade de detec��o
% Pfa - Probabilidade de falso alarme
% T_H0 - Vetor das estat�sticas sob H0 (canal ocioso)
% T_H1 - Vetor das estat�sticas sob H1 (canal ocupado)
%%

function [Pd,Pfa,T_H0,T_H1] = func_calc_Pd_Pfa_Phase_Dif_vf(Nmc,Nsimb,Channel,SNRdb,M,Sepec_Sens_Technique,Thre,nbins,Sample_factor,Noise_uncertainty_db,Error_fc,parametro,ParametroName,CannelParam,INParam)

%%
% INICIALIZA VETORES PARA ARMAZENAR ESTAT�STICAS SOB H0 (CANAL OCIOSO) E H1 (CANAL OCUPADO)

T_H0 = zeros(1,Nmc/2);
T_H1 = zeros(1,Nmc/2);
%%

%%
% CONFIGURA��ES DO SISTEMA

% Eventos de Montecarlo
% Nmc = 1000;

% N�mero de s�mbolos transmitidos
% Nsimb = 1000;

% Canal Rayleigh lento (1) ou gaussino (0)
% Channel = 1;

% SNR na recep��o
% SNRdb = -20;

% Ordem da modula��o MPSK considerada
% M = 2;
%%

%%
% DEFINI��ES DAS TAXAS ENVOLVIDAS

% Taxa de bit
Rb = 1;

% Taxa de s�mbolo
Rs = Rb ./ log2(M);

% Tempo de s�mbolo
Tsimb = 1/Rs;

% Fator de amostragem (quantas vezes a frequ�ncia de amostragem � maior que a taxa de s�mbolo)
K = Sample_factor;

% Taxa de amostragem
Fs = K*Rs;

% Tempo de amostragem
Ts = 1/Fs;

% Frequ�ncia do oscilador local
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
% INCERTEZA DE RU�DO EM ESCALA LINEAR

Noise_uncertainty = 10.^(Noise_uncertainty_db/10);
%%

%%
% LOOP MONTE CARLO

for i = 1 : Nmc
    
    %%
    % DEFINE SE � HIP�TESE H0 (CANAL OCIOSO) OU H1 (CANAL OCUPADO)
    
    if i <= Nmc/2
        TxEnable = 0;
    else
        TxEnable = 1;
    end
    %%
    
    %%
    % GERA��O DOS S�MBOLOS TRANSMITIDOS
    
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
    % GERA��O DOS DADOS AMOSTRADOS
    
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
    % CRIA��O DO VETOR DE TEMPO PARA O SINAL MODULADO
    
    t = [0: 1/Fs : (length(Data_Kx) - 1)*1/Fs];
    % cosseno = sqrt(2/Tsimb)*cos(2*pi*fosc.*t);
    % seno = sqrt(2/Tsimb)*sin(2*pi*fosc.*t);
    %%
    
    %%
    % MODULA��O COMPLEXA COM O OSCILADOR LOCAL (TRANSLA��O EM FREQU�NCIA)
    
    % Data_Mod = real(Data_Kx) .*cosseno + 1i*imag(Data_Kx).*seno;
    
    % S�mbolo multiplicado pela exponencial complexa e^(j2*pi*fc/fst)
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
    % NORMALIZA��O DE ENERGIA E ADI��O DE RU�DO AWGN COM INCERTEZA DE RU�DO
    
    Pot_sinal = sum(real(Data_Mod).^2 + imag(Data_Mod).^2)./length(Data_Mod);
    Data_Mod = Data_Mod/sqrt(Pot_sinal)*sqrt(10.^(SNRdb/10));
    Noise = sqrt(unifrnd(1/Noise_uncertainty,Noise_uncertainty))*(1/sqrt(2)*randn(1,length(Data_Mod)) + 1i* 1/sqrt(2)*randn(1,length(Data_Mod)));
    % Data_rx = TxEnable*Data_Mod + Noise;
    
    % Aplica��o do canal (fading) e ru�do impulsivo se necess�rio
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
    % ESTIMA��O DAS FASES INSTANT�NEAS
    
    %fi_n = mod(phase(Data_rx),2*pi); % Estava demorando muito
    %fi_n = mod(atan2(imag(Data_rx),real(Data_rx)),2*pi);
    fi_n = atan2(imag(Data_rx),real(Data_rx));
    %%
    
    %%
    % CALCULA A DIFEREN�A DE FASE ENTRE AMOSTRAS CONSECUTIVAS
    
    theta_n = mod((fi_n(2:end)-fi_n(1:end-1)),2*pi);
    
    %Plota histogramas da fase do sinal recebido e da diferen�a de fase do sinal recebido
    % figure; hist(fi_n,50);legend('pdf de phi_n');
    % figure; hist(theta_n,50);legend('pdf de theta');
    %%
    
    %%
    % CALCULA A PDF EMP�RICA DA DIFEREN�A DE FASE
    
    [xpdf_empirical,ypdf_empirical] = fcn_coder_pdf(theta_n, nbins);
    % figure; stem(xpdf_empirical,ypdf_empirical,'--*b')
    % var(ypdf_empirical);
    %%
    
    %%
    % PLOTANDO AS PDFS TE�RICAS
    
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
    % xlabel('Diferen�a de fase [rad]')
    % ylabel('Pr')
    %%
    
    %%
    % C�LCULO DA VARI�VEL DE DECIS�O COM BASE NA T�CNICA DE DETEC��O ESCOLHIDA
    
    [fmax, pxmax]= max(ypdf_empirical);
    
    % Detec��o de energia
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
        
        % Fallback: detec��o com pot�ncia fracion�ria
    else
        T = sum(abs(Data_rx).^1.1);
    end
    %%
    
    %%
    % ENCONTRANDO A THRE
    
    % Thre = 0.1127/sqrt(N)*(qfuncinv(Pfa))
    %%
    
    %%
    % ARMAZENA ESTAT�STICA E CONTA PD OU PFA
    
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
