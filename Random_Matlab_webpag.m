%%
% GERA NÚMEROS ALEATÓRIOS A PARTIR DE UMA DISTRIBUIÇÃO DE PROBABILIDADE DEFINIDA PELO USUÁRIO

% Função: Random_Matlab_Webpag

% Sintaxe:
% x = randpdf(p, px, dim)
% randpdf(p, px, dim)

% Descrição:
% x = randpdf(p, px, dim) retorna uma matriz de números aleatórios a partir
% da distribuição de densidade de probabilidade definida em p e px.
% p são os valores da densidade (eixo y) e px os valores correspondentes
% da pdf (eixo x). p e px devem ter o mesmo comprimento.
% dim define as dimensões da matriz de saída. Por exemplo, dim = [100 3] define
% uma matriz 100x3 (com 300 números aleatórios).

% Atenção: Isto não é um gerador de números realmente aleatórios,
% mas apenas uma transformação de números pseudorrandômicos uniformemente
% distribuídos para a pdf desejada!

% Exemplo 1:
% Geração de números aleatórios distribuídos normalmente.
% Esta não é uma distribuição normal típica, pois está limitada
% dos dois lados, ou seja, 0 < px < 80.

% px = 0:80;
% p = 1./(10*sqrt(2*pi))*exp((-(px-40).^2)./(2*10^2));
% randpdf(p,px,[10000,1])

% Exemplo 2:
% Geração usando uma pdf definida pelo usuário.

% px = [1 2 3 4 5 6 7 8 9];
% p = [0 1 3 0 0 4 5 4 0];
% randpdf(p,px,[50000,1])

% Entradas:
% p - densidade de probabilidade,
% px - valores correspondentes à densidade de probabilidade,
% dim - dimensão da matriz de saída.

% SAÍDA:
% x - números aleatórios. Execute a função sem saída para visualizar gráficos.

% Por Adam Niesony, Universidade de Tecnologia de Opole, Polônia
%%

function x = Random_Matlab_webpag(p,px,dim)

%%
% VERIFICA O NÚMERO DE ARGUMENTOS PASSADOS

error(nargchk(3, 3, nargin))
%%

%%
% ORGANIZA EM VETORES COLUNA E CALCULA A PDF

px = px(:);

% Normaliza para integral = 1 (área sob curva)
p = p(:)./trapz(px,p(:));
%%

%%
% REALIZA UMA INTERPOLAÇÃO PARA TER UMA REPRESENTAÇÃO MAIS SUAVE DA PDF

% Cria um vetor x interpolado com 10.000 pontos igualmente espaçados
pxi = [linspace(min(px),max(px),10000)]';

% Interpola os valores da pdf para os pontos criados
pi = interp1(px,p,pxi,'linear');
%%

%%
% CALCULA A FUNÇÃO DISTRIBUIÇÃO ACUMULADA (CDF)

% A CDF é a integral acumulada da pdf
cdfp = cumtrapz(pxi,pi);
%%

%%
% REMOVE PARTES CONSTANTES NA CDF PARA EVITAR PROBLEMAS NA INVERSÃO (PONTOS ONDE A CDF É "PLANA")

ind = [true; not(diff(cdfp)==0)];
cdfp = cdfp(ind);
pi = pi(ind);
pxi = pxi(ind);
%%

%%
% GERA NÚMEROS UNIFORMEMENTE DISTRIBUÍDOS ENTRE 0 E 1

uniformDistNum = rand(dim);
%%

%%
% USA A CDF PARA TRANSFORMAR NÚMEROS UNIFORMES EM NÚMEROS COM A DISTRIBUIÇÃO DESEJADA

% Técnica da inversão da CDF
userDistNum = interp1(cdfp,pxi,uniformDistNum(:)','linear');
%%

%%
% PLOTAGEM DOS GRÁFICOS

% Se não houver argumento de saída, plota gráficos ilustrando o resultado,
% caso tenha argumento de saída, retorna matriz na dimensão pedida

if nargout == 0
    subplot(3,4,[1 2 5 6])
    
    % Histograma
    [n,xout] = hist(userDistNum,50);
    
    % Normaliza o Histograma
    n = n./sum(n)./(xout(2)-xout(1));
    bar(xout,n)
    hold on
    
    % Plota a PDF original
    plot(pxi, pi./trapz(pxi,pi),'r')
    hold off
    legend('pdf dos números gerados', 'pdf original')
    
    subplot(3,4,[3 4 7 8])
    plot(pxi, cdfp,'g')
    ylim([0 1])
    legend('cdf da pdf de entrada')
    
    subplot(3,4,[9:12])
    plot(userDistNum)
    legend('números gerados')
else
    x = reshape(userDistNum,dim);
end
%%