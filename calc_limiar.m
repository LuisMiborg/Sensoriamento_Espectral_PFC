%%
% CALCULA UM VETOR DE LIMIARES (THRESHOLDS) EM ESCALA LOGAR�TMICA

% Fun��o: calc_limiar

% Autores: Andr� Ant�nio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberl�ndia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descri��o:
% Gera um vetor de limiares distribu�dos de forma logar�tmica
% entre 'limiar_min' e 'limiar_max'. � �til para varrer diferentes valores
% de limiar, por exemplo em experimentos de detec��o, varrendo thresholds
% para construir curvas ROC ou similares.

% Entradas:
% N_pts - N�mero de pontos (quantidade de limiares a serem calculados)
% limiar_min - Valor m�nimo do limiar
% limiar_max - Valor m�ximo do limiar

% Sa�das:
% limiar - Vetor (1 x N_pts) contendo os valores dos limiares calculados

function [limiar] = calc_limiar(N_pts,limiar_min, limiar_max)

% Cada limiar � calculado em escala logar�tmica
% usando log(i) / log(N_pts) para obter uma distribui��o entre 0 e 1,
% multiplicando pela faixa (limiar_max - limiar_min) e somando o m�nimo
for i = 1 : N_pts
    
    limiar(i) = log(i) / log(N_pts) * (limiar_max - limiar_min) + limiar_min;
    
end