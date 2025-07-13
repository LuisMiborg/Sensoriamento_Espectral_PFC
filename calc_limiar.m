%%
% CALCULA UM VETOR DE LIMIARES (THRESHOLDS) EM ESCALA LOGARÍTMICA

% Função: calc_limiar

% Autores: André Antônio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberlândia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descrição:
% Gera um vetor de limiares distribuídos de forma logarítmica
% entre 'limiar_min' e 'limiar_max'. É útil para varrer diferentes valores
% de limiar, por exemplo em experimentos de detecção, varrendo thresholds
% para construir curvas ROC ou similares.

% Entradas:
% N_pts - Número de pontos (quantidade de limiares a serem calculados)
% limiar_min - Valor mínimo do limiar
% limiar_max - Valor máximo do limiar

% Saídas:
% limiar - Vetor (1 x N_pts) contendo os valores dos limiares calculados

function [limiar] = calc_limiar(N_pts,limiar_min, limiar_max)

% Cada limiar é calculado em escala logarítmica
% usando log(i) / log(N_pts) para obter uma distribuição entre 0 e 1,
% multiplicando pela faixa (limiar_max - limiar_min) e somando o mínimo
for i = 1 : N_pts
    
    limiar(i) = log(i) / log(N_pts) * (limiar_max - limiar_min) + limiar_min;
    
end