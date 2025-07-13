%%
% GERA AMOSTRAS ALEATÓRIAS SEGUNDO DIFERENTES MODELOS DE DESVANECIMENTO (FADING)

% Função: Random_Samples_Fading_model_vf

% Autores: André Antônio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberlândia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descrição:
% Gera amostras aleatórias segundo diferentes modelos estatísticos
% de desvanecimento (fading) de canal, tais como Rayleigh, Rice, Nakagami-m,
% Alpha-Mu, Kappa-Mu e Eta-Mu. Calcula e retorna a PDF teórica
% do modelo escolhido e o vetor de amplitudes normalizadas usado para gerar
% as amostras.

% Entradas:
% fading_type - Tipo de modelo de desvanecimento:
% 0 = Rayleigh
% 1 = Rice
% 2 = Nakagami-m
% 3 = Alpha-Mu
% 4 = Kappa-Mu
% 5 = Eta-Mu
% nsamples - Número de amostras aleatórias a serem geradas.
% parametro - Valor do parâmetro principal do modelo (ex.: m para Nakagami)
% ParametroName - String indicando qual parâmetro está sendo variado,
% ou 'default' para usar valores padrão definidos no struct.
% CannelParam  - Struct contendo os parâmetros default de cada modelo.

% Saídas:
% x - Vetor (1 x nsamples) com as amostras aleatórias geradas segundo
% o modelo de desvanecimento especificado.
% r - Vetor com os valores de amplitude normalizada usados para calcular a
% PDF (ex.: de 0 a 10, com passo 0.01).
% fr - Vetor com a PDF teórica do modelo de desvanecimento avaliada nos
% pontos de r.
%%

function [x,r,fr] = Random_Samples_Fading_model_vf(fading_type,nsamples,parametro,ParametroName,CannelParam)

%% 
% PDF RAYLEIGH

p00 = 0.707;
fr_Rayleigh = @(p,r) r/(p(1).^2) .*exp(-(r.^2/(2*p(1).^2)));
% p(1) = sig ... rc = sqrt(2*sig.^2)
p00 = 1;
fr_Rayleigh = @(p,r) 1/p(1)*2*(r/p(1)) .*exp(-(r/p(1)).^2);
% p(1) = rc
%%

%% 
% PDF RICE

p01 = [0.707 0.000000001];
fr_Rice = @(p,r) r/(p(1).^2) .*exp(-((r.^2 + p(2).^2)/(2*p(1).^2))).*besseli(0,(r*p(2))/p(1).^2);
% p(1) = sig ; p(2) = v ... k = v.^2/(2*sig.^2)

if strcmp(ParametroName,'default')
    p01 = [CannelParam.Rice.kappa 1];
else
    p01 = [parametro 1];
end

fr_Rice = @(p,r) 1/p(2)*2*(p(1)+1)*(r/p(2)).*exp(-p(1)-(p(1)+1)*(r/p(2)).^2).*besseli(0,2*sqrt(p(1)*(1+p(1)))*(r/p(2)));
% p(1) = k ; p(2) = rc
%%

%% 
% PDF NAKAGAMI-M

if strcmp(ParametroName,'default')
    p02 = [CannelParam.Naka.mu  1];
else
    p02 = [parametro 1];
end

fr_Nakagami = @(p,r)( 1/p(2)*2*p(1).^p(1) * (r/p(2)).^(2*p(1)-1) .* exp(-p(1)*(r/p(2)).^2))/gamma(p(1));
% p(1) = m; p(2)= rc
%%

%% 
% PDF ALPHAMU

if strcmp(ParametroName,'default')
    p03 = [CannelParam.AlphaMu.alpha CannelParam.AlphaMu.mu 1];
elseif strcmp(ParametroName,'alpha')
    p03 = [parametro CannelParam.AlphaMu.mu 1];
elseif strcmp(ParametroName,'mu')
    p03 = [CannelParam.AlphaMu.alpha parametro 1];
end

% p03 = [parametro 2 1];
fr_Alphamu = @(p,r)(p(1) * p(2).^p(2) * r.^(p(1)*p(2)-1))/(p(3).^(p(1)*p(2))*gamma(p(2))).*exp(-p(2)*r.^p(1)/(p(3).^p(1)));
%p(1)= alpha; p(2)= mu; p(3)= rc
%%

%% 
% PDF KAPPA-MU

% p04 =[0.00001 2 1];
% p04 =[10 2 1];
% p04 = [parametro 2 1]; %Variando kappa ROC 2
% p04 = [parametro 1 1]; %Variando kappa ROC 2
% p04 = [0.00001 parametro 1]; %Variando mu ROC 3

if strcmp(ParametroName,'default')
    p04 =[CannelParam.KappaMu.kappa CannelParam.KappaMu.mu 1];
elseif strcmp(ParametroName,'kappa')
    p04 =[parametro CannelParam.KappaMu.mu 1];
elseif strcmp(ParametroName,'mu')
    p04 =[CannelParam.KappaMu.kappa parametro 1];
end

fr_Kappamu = @(p,r)(1/p(3)*2*p(2)*(1+p(1)).^((p(2)+1)/2))/(p(1).^((p(2)-1)/2)*exp(p(2)*p(1))).*(r/p(3)).^(p(2)).*exp(-p(2).*(1+p(1)).*(r/p(3)).^2).*besseli((p(2)-1),(2*p(2)*sqrt(p(1)*(1+p(1))).*(r/p(3))));
%p(1) = kappa; p(2) = mu; p(3) = rc
%%

%% 
% PDF ETA-MU

if strcmp(ParametroName,'default')
    p05 =[CannelParam.EtaMu.eta CannelParam.EtaMu.mu 1];
elseif strcmp(ParametroName,'eta')
    p05 =[parametro CannelParam.EtaMu.mu 1];
elseif strcmp(ParametroName,'mu')
    p05 =[CannelParam.EtaMu.eta parametro 1];
end

% p03 = [parametro 1/2 1];
fr_EtaMu = @(p,r) abs((1/p(3))*((4*sqrt(pi)*(p(2)^(p(2)+0.5))*((2+p(1)^(-1)+p(1))/4)^p(2))/(gamma(p(2))*((p(1)^(-1)-p(1))/4)^(p(2)-0.5))).*((r/p(3)).^(2*p(2))).*exp(-2*p(2)*((2+p(1)^(-1)+p(1))/4).*(r/p(3)).^2).*besseli((p(2)-0.5),(2*p(2)*((p(1)^(-1)-p(1))/4).*(r/p(3)).^2)));

% Expressao da EtaMu utilizando simplificacao da AlphaEtaKappaMu
% pdfEnvoltoriaPhase = @(theta,r,alfa,eta,kappa,mu,p,q,rc) ((alfa*(mu^2)*p*(eta+1)^2*(kappa+1)^((mu/2)+1).*r.^((alfa/2)*(mu+2)-1)*(abs(sin(theta))).^(mu/(p+1)).*(abs(cos(theta))).^(mu*p/(p+1)))/...
% (2*eta*((p+1)^2)*((kappa/(eta*q+1))^((mu/2)-1))*(eta*q)^((mu*p)/(2*(p+1))-0.5)*(rc^((alfa/2)*(mu+2)))*exp((kappa*mu*(eta+1)*(q*p+1))/((p+1)*(eta*q+1))))).*...
% exp(-((mu*(eta+1)*(kappa+1)*((eta*(sin(theta)).^2)+(p*(cos(theta)).^2)))/(eta*(p+1))).*((r/rc).^alfa)).*...
% exp(((2*mu*(eta+1)*cos(theta-atan((1/p)*((eta/q))^(0.5))))/(eta*(p+1)))*sqrt((eta*kappa*(kappa+1)*(eta+q*(p^2)))/(eta*q+1)).*(r/rc).^(alfa/2)).*...
% (besseli(((mu/(p+1))-1),(((2*mu*(eta+1)*abs(sin(theta)))/(p+1))*sqrt((kappa*(kappa+1))/(eta*q+1)).*(r/rc).^(alfa/2))).*besseli((((mu*p)/(p+1))-1),(((2*mu*p*(eta+1)*abs(cos(theta)))/(eta*(p+1)))*sqrt((eta*kappa*q*(kappa+1))/(eta*q+1)).*(r/rc).^(alfa/2)))./...
% (cosh(((2*mu*(eta+1)*sin(theta))/(p+1))*sqrt((kappa*(kappa+1))/(eta*q+1)).*(r/rc).^(alfa/2)).*cosh(((2*mu*p*(eta+1)*cos(theta))/(eta*(p+1)))*sqrt((eta*kappa*q*(kappa+1))/(eta*q+1)).*(r/rc).^(alfa/2))));
% fr_EtaMu = @(p,r) (integral(@(theta) pdfEnvoltoriaPhase(theta,r,2,p(1),0.00001,p(2),1,1,p(3)),0,2*pi,'ArrayValued',true));
%%

%% 
% PDF ALPHAETAKAPPAMU

% p06 = [2 1 0.0001 1 1 1 1];
% pdfEnvoltoriaPhase = @(theta,r,alfa,eta,kappa,mu,p,q,rc) ((alfa*(mu^2)*p*(eta+1)^2*(kappa+1)^((mu/2)+1).*r.^((alfa/2)*(mu+2)-1)*(abs(sin(theta))).^(mu/(p+1)).*(abs(cos(theta))).^(mu*p/(p+1)))/...
% (2*eta*((p+1)^2)*((kappa/(eta*q+1))^((mu/2)-1))*(eta*q)^((mu*p)/(2*(p+1))-0.5)*(rc^((alfa/2)*(mu+2)))*exp((kappa*mu*(eta+1)*(q*p+1))/((p+1)*(eta*q+1))))).*...
% exp(-((mu*(eta+1)*(kappa+1)*((eta*(sin(theta)).^2)+(p*(cos(theta)).^2)))/(eta*(p+1))).*((r/rc).^alfa)).*...
% exp(((2*mu*(eta+1)*cos(theta-atan((1/p)*((eta/q))^(0.5))))/(eta*(p+1)))*sqrt((eta*kappa*(kappa+1)*(eta+q*(p^2)))/(eta*q+1)).*(r/rc).^(alfa/2)).*...
% (besseli(((mu/(p+1))-1),(((2*mu*(eta+1)*abs(sin(theta)))/(p+1))*sqrt((kappa*(kappa+1))/(eta*q+1)).*(r/rc).^(alfa/2))).*besseli((((mu*p)/(p+1))-1),(((2*mu*p*(eta+1)*abs(cos(theta)))/(eta*(p+1)))*sqrt((eta*kappa*q*(kappa+1))/(eta*q+1)).*(r/rc).^(alfa/2)))./...
% (cosh(((2*mu*(eta+1)*sin(theta))/(p+1))*sqrt((kappa*(kappa+1))/(eta*q+1)).*(r/rc).^(alfa/2)).*cosh(((2*mu*p*(eta+1)*cos(theta))/(eta*(p+1)))*sqrt((eta*kappa*q*(kappa+1))/(eta*q+1)).*(r/rc).^(alfa/2))));
% fr_AlphaEtaKappaMu = @(p,r) (integral(@(theta) pdfEnvoltoriaPhase(theta,r,p(1),p(2),p(3),p(4),p(5),p(6),p(7)),0,2*pi,'ArrayValued',true));
% % p(1) = alpha; p(2) = eta; p(3) = kappa; p(4) = mu; p(5) = p; p(6) = q; % p(7) = rc;
%%

%% 
% DEFINE VETOR DE VALORES PARA AVALIAR A PDF

r = 0:0.01:10;
%%

%% 
% SELECIONA MODELO DE FADING CONFORME TIPO ESPECIFICADO

if fading_type == 0
    fr = fr_Rayleigh(p00,r);
elseif fading_type == 1
    fr = fr_Rice(p01,r);
elseif fading_type == 2
    fr = fr_Nakagami(p02,r);
elseif fading_type == 3
    fr= fr_Alphamu(p03,r);
elseif fading_type == 4
    fr= fr_Kappamu(p04,r);
elseif fading_type ==5
    fr= fr_EtaMu(p05,r);
end

% figure;
% plot(r,fr,'k')
% hold on;
% xlabel('Normalized Envelope Level')
% ylabel('PDF')
% plot(r,fr_Rayleigh_lin,'k');
% plot(r,fr_AlphaMu_lin,'--y')
% plot(r,fr_KappaMu_lin,'*b')
%%

%% 
% GERA AMOSTRAS ALEATÓRIAS A PARTIR DA PDF ESCOLHIDA
x = Random_Matlab_webpag(abs(fr),r,[1,nsamples]);
%%

% figure;
% hist(x);

