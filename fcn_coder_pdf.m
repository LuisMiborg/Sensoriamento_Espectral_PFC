%%
% CALCULA UMA ESTIMATIVA DA PDF A PARTIR DE DADOS AMOSTRAIS

% Função: fcn_coder_pdf

% Autores: André Antônio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberlândia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descrição:
% Estima a PDF de um conjunto de dados amostrais, utilizando a
% técnica de histograma. O histograma é normalizado para que a área sob a
% curva seja igual a 1, aproximando assim a PDF dos dados.

% Entradas:
% data - Vetor com os dados amostrais
% bins - Número de bins (caixas) para o histograma

% Saídas:
% xo - Centros dos bins (eixo x)
% yo - Valores normalizados do histograma (pdf estimada no eixo y)
%%

function [xo,yo] = fcn_coder_pdf(data, bins)

% Calcula o histograma dos dados
[n,xout] = hist(data,bins);

% Normaliza o histograma para ter área unitária (aproxima uma PDF)
% dividindo pelo número total de amostras e pelo tamanho do bin
n = n./sum(n)./(xout(2)-xout(1));

% Salva os centros dos bins no eixo x
xo = xout;

% Salva a densidade normalizada no eixo y
yo = n;

% Calcula a área sob a curva (para verificar se é ~1)
Area = sum(n).*(xout(2)-xout(1)); % esse valor deve ser próximo de 1 se a normalização estiver correta