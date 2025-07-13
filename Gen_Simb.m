%%
% FUN��O PARA GERAR UM VETOR DE N S�MBOLOS MODULADOS, CONFORME UMA CONSTELA��O DE ORDEM M

% Fun��o: Gen_Simb

% Autores: Andr� Ant�nio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberl�ndia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descri��o:
% Gera um vetor de N s�mbolos modulados a partir de uma
% constela��o de ordem M, normalizada em energia unit�ria.

% Entradas:
% M - Ordem da constela��o (ex: 2 para BPSK, 4 para QPSK, etc.)
% N - N�mero de s�mbolos a serem gerados

% Sa�das:
% Simbolos - Vetor 1 x N contendo os s�mbolos modulados, sorteados
% aleatoriamente da constela��o normalizada
%%

function [Simbolos] = Gen_Simb(M,N)

%% 
%DEFINI��O DA CONSTELA��O

% Caso BPSK com s�mbolos +1 e -1
% aux = modem.qammod('M',M);
Vec_Simbol = [1 -1];

% Recebe o vetor da constela��o MQAM
Simbolos_Mqam = Vec_Simbol;

% Calcula a energia m�dia dos s�mbolos da constela��o para posterior normaliza��o
Energy = (sum(real(Simbolos_Mqam).^2) + sum(imag(Simbolos_Mqam).^2))/ length(Simbolos_Mqam);

% Normaliza o vetor de s�mbolos para energia unit�ria, ou seja, E[|S|^2] = 1
cte_norm = sqrt(Energy);
Simbolos_Mqam = Simbolos_Mqam / sqrt(Energy);

% Gera um vetor de N s�mbolos aleat�rios, escolhidos da constela��o MQAM normalizada
SOFDM = Simbolos_Mqam(randi([1,length(Simbolos_Mqam)],1,N));

% Retorna o vetor de s�mbolos gerado
Simbolos = SOFDM;
%%
