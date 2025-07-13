%%
% GERA N�MEROS ALEAT�RIOS A PARTIR DE UMA DISTRIBUI��O DE PROBABILIDADE DEFINIDA PELO USU�RIO

% Fun��o: Random_Matlab_Webpag

% Sintaxe:
% x = randpdf(p, px, dim)
% randpdf(p, px, dim)

% Descri��o:
% x = randpdf(p, px, dim) retorna uma matriz de n�meros aleat�rios a partir
% da distribui��o de densidade de probabilidade definida em p e px.
% p s�o os valores da densidade (eixo y) e px os valores correspondentes
% da pdf (eixo x). p e px devem ter o mesmo comprimento.
% dim define as dimens�es da matriz de sa�da. Por exemplo, dim = [100 3] define
% uma matriz 100x3 (com 300 n�meros aleat�rios).

% Aten��o: Isto n�o � um gerador de n�meros realmente aleat�rios,
% mas apenas uma transforma��o de n�meros pseudorrand�micos uniformemente
% distribu�dos para a pdf desejada!

% Exemplo 1:
% Gera��o de n�meros aleat�rios distribu�dos normalmente.
% Esta n�o � uma distribui��o normal t�pica, pois est� limitada
% dos dois lados, ou seja, 0 < px < 80.

% px = 0:80;
% p = 1./(10*sqrt(2*pi))*exp((-(px-40).^2)./(2*10^2));
% randpdf(p,px,[10000,1])

% Exemplo 2:
% Gera��o usando uma pdf definida pelo usu�rio.

% px = [1 2 3 4 5 6 7 8 9];
% p = [0 1 3 0 0 4 5 4 0];
% randpdf(p,px,[50000,1])

% Entradas:
% p - densidade de probabilidade,
% px - valores correspondentes � densidade de probabilidade,
% dim - dimens�o da matriz de sa�da.

% SA�DA:
% x - n�meros aleat�rios. Execute a fun��o sem sa�da para visualizar gr�ficos.

% Por Adam Niesony, Universidade de Tecnologia de Opole, Pol�nia
%%

function x = Random_Matlab_webpag(p,px,dim)

%%
% VERIFICA O N�MERO DE ARGUMENTOS PASSADOS

error(nargchk(3, 3, nargin))
%%

%%
% ORGANIZA EM VETORES COLUNA E CALCULA A PDF

px = px(:);

% Normaliza para integral = 1 (�rea sob curva)
p = p(:)./trapz(px,p(:));
%%

%%
% REALIZA UMA INTERPOLA��O PARA TER UMA REPRESENTA��O MAIS SUAVE DA PDF

% Cria um vetor x interpolado com 10.000 pontos igualmente espa�ados
pxi = [linspace(min(px),max(px),10000)]';

% Interpola os valores da pdf para os pontos criados
pi = interp1(px,p,pxi,'linear');
%%

%%
% CALCULA A FUN��O DISTRIBUI��O ACUMULADA (CDF)

% A CDF � a integral acumulada da pdf
cdfp = cumtrapz(pxi,pi);
%%

%%
% REMOVE PARTES CONSTANTES NA CDF PARA EVITAR PROBLEMAS NA INVERS�O (PONTOS ONDE A CDF � "PLANA")

ind = [true; not(diff(cdfp)==0)];
cdfp = cdfp(ind);
pi = pi(ind);
pxi = pxi(ind);
%%

%%
% GERA N�MEROS UNIFORMEMENTE DISTRIBU�DOS ENTRE 0 E 1

uniformDistNum = rand(dim);
%%

%%
% USA A CDF PARA TRANSFORMAR N�MEROS UNIFORMES EM N�MEROS COM A DISTRIBUI��O DESEJADA

% T�cnica da invers�o da CDF
userDistNum = interp1(cdfp,pxi,uniformDistNum(:)','linear');
%%

%%
% PLOTAGEM DOS GR�FICOS

% Se n�o houver argumento de sa�da, plota gr�ficos ilustrando o resultado,
% caso tenha argumento de sa�da, retorna matriz na dimens�o pedida

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
    legend('pdf dos n�meros gerados', 'pdf original')
    
    subplot(3,4,[3 4 7 8])
    plot(pxi, cdfp,'g')
    ylim([0 1])
    legend('cdf da pdf de entrada')
    
    subplot(3,4,[9:12])
    plot(userDistNum)
    legend('n�meros gerados')
else
    x = reshape(userDistNum,dim);
end
%%