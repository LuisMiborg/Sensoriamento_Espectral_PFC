%%
% CALCULA UMA ESTIMATIVA DA PDF A PARTIR DE DADOS AMOSTRAIS

% Fun��o: fcn_coder_pdf

% Autores: Andr� Ant�nio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberl�ndia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descri��o:
% Estima a PDF de um conjunto de dados amostrais, utilizando a
% t�cnica de histograma. O histograma � normalizado para que a �rea sob a
% curva seja igual a 1, aproximando assim a PDF dos dados.

% Entradas:
% data - Vetor com os dados amostrais
% bins - N�mero de bins (caixas) para o histograma

% Sa�das:
% xo - Centros dos bins (eixo x)
% yo - Valores normalizados do histograma (pdf estimada no eixo y)
%%

function [xo,yo] = fcn_coder_pdf(data, bins)

% Calcula o histograma dos dados
[n,xout] = hist(data,bins);

% Normaliza o histograma para ter �rea unit�ria (aproxima uma PDF)
% dividindo pelo n�mero total de amostras e pelo tamanho do bin
n = n./sum(n)./(xout(2)-xout(1));

% Salva os centros dos bins no eixo x
xo = xout;

% Salva a densidade normalizada no eixo y
yo = n;

% Calcula a �rea sob a curva (para verificar se � ~1)
Area = sum(n).*(xout(2)-xout(1)); % esse valor deve ser pr�ximo de 1 se a normaliza��o estiver correta