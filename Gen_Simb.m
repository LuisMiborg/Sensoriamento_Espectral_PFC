%%
% FUNÇÃO PARA GERAR UM VETOR DE N SÍMBOLOS MODULADOS, CONFORME UMA CONSTELAÇÃO DE ORDEM M

% Função: Gen_Simb

% Autores: André Antônio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberlândia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descrição:
% Gera um vetor de N símbolos modulados a partir de uma
% constelação de ordem M, normalizada em energia unitária.

% Entradas:
% M - Ordem da constelação (ex: 2 para BPSK, 4 para QPSK, etc.)
% N - Número de símbolos a serem gerados

% Saídas:
% Simbolos - Vetor 1 x N contendo os símbolos modulados, sorteados
% aleatoriamente da constelação normalizada
%%

function [Simbolos] = Gen_Simb(M,N)

%% 
%DEFINIÇÃO DA CONSTELAÇÃO

% Caso BPSK com símbolos +1 e -1
% aux = modem.qammod('M',M);
Vec_Simbol = [1 -1];

% Recebe o vetor da constelação MQAM
Simbolos_Mqam = Vec_Simbol;

% Calcula a energia média dos símbolos da constelação para posterior normalização
Energy = (sum(real(Simbolos_Mqam).^2) + sum(imag(Simbolos_Mqam).^2))/ length(Simbolos_Mqam);

% Normaliza o vetor de símbolos para energia unitária, ou seja, E[|S|^2] = 1
cte_norm = sqrt(Energy);
Simbolos_Mqam = Simbolos_Mqam / sqrt(Energy);

% Gera um vetor de N símbolos aleatórios, escolhidos da constelação MQAM normalizada
SOFDM = Simbolos_Mqam(randi([1,length(Simbolos_Mqam)],1,N));

% Retorna o vetor de símbolos gerado
Simbolos = SOFDM;
%%
