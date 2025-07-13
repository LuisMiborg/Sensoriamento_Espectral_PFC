%%
% SIMULA O EFEITO DE RU�DO IMPULSIVO NAS AMOSTRAS COLETADAS

% Fun��o: Gen_Impulsive_Noise

% Autores: Andr� Ant�nio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberl�ndia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descri��o:
% Simula o efeito do ru�do impulsivo sobre sinais complexos coletados
% em m receptores (ou linhas), cada um com n amostras.
% Permite gerar dois tipos de ru�do impulsivo (Op_IN_Tipe = 1 ou 2).

% Entradas:
% m - N�mero de linhas (ex: sensores ou receptores)
% n - N�mero de amostras por linha
% Nb - N�mero de bursts (pulsos) de ru�do impulsivo (tipo 1)
% Ns - Largura de cada burst de ru�do impulsivo (tipo 1)
% K - Rela��o entre pot�ncia do ru�do impulsivo e do ru�do gaussiano
% Pin - Probabilidade do ru�do impulsivo estar presente
% Pcr - Probabilidade de um receptor ser afetado pelo ru�do impulsivo
% A, B - Par�metros do ru�do lognormal (tipo 2)
% Beta - Par�metro da distribui��o exponencial para espa�amento (tipo 2)
% Op_IN_Tipe Tipo do ru�do impulsivo:
% 1: bursts localizados (pulsos)
% 2: rajadas com espa�amento exponencial
% SNR - Rela��o sinal-ru�do em dB (usado para calcular sigma_i no tipo 1)

% Sa�das:
% VI - Matriz m x n de ru�do impulsivo complexo
%%

function [VI] = Gen_Impulsive_Noise(m,n,Nb,Ns,K,Pin,Pcr,A,B,Beta,Op_IN_Tipe,SNR)

%%
% RU�DO IMPULSIVO TIPO 1

if Op_IN_Tipe == 1
    
    % Calcula o desvio padr�o do ru�do impulsivo baseado na rela��o SNR e fator K
    Sigma_i = sqrt( K / (10^(SNR/10)));
    
    % Inicializa matriz Z com zeros para indicar posi��es dos pulsos
    Z = zeros(m,n);
    
    % Para cada linha (para cada r�dio ou sensor)
    for u = 1 : m
        
        % Inicializa vetor de zeros para marcar posi��es de pulso
        g = zeros(1,n);
        
        % Defini��o dos intervalos entre os pulsos de ru�do impulsivo
        t(2 * Nb) = 0;
        
        % Primeiro intervalo aleat�rio para posicionar o primeiro pulso
        t(1) = randi([0 (floor(n /Nb))]);
        t(2) = t(1) + Ns - 1; % Fim do primeiro pulso
        
        % Se houver m�ltiplos pulsos (Nb > 1), define os intervalos subsequentes
        if Nb > 1
            for r = 1 : Nb -1
                t(2*r + 1) = t(2*r) + randi([0 (floor(n /Nb))]) + Ns;
                t(2*r + 2) = t(2*r + 1) + Ns - 1;
            end
        end
        
        % Copia os intervalos para z
        z = t;
        
        % Marca as posi��es dos pulsos no vetor g (com 1)
        for r = 0 : Nb -1
            for w = z(2*r+1) : z(2*r +2)
                g(mod(w,n)+1)= 1;  % Usa m�dulo para circular no vetor
            end
        end
        
        % Preenche a linha da matriz Z com os pulsos marcados
        Z(u,:) = g;
    end
    
    % Determina quais r�dios ser�o afetados pelo ru�do impulsivo, via vari�vel binomial
    for i = 1 : m
        aux_bino(i) = binornd(1,Pcr);
    end
    
    % Remove ru�do dos r�dios n�o afetados
    G = diag(aux_bino) * Z;
    
    N = zeros(m,n);
    
    % Gera ru�do Gaussiano complexo branco
    N = randn(m,n) +  1i * randn(m,n);
    
    % Aplica ru�do somente nas posi��es indicadas em G
    V2 = N.*G;
    
    % Normaliza e aplica fator Sigma_i para cada r�dio
    for i = 1 : m
        
        % Calcula energia do ru�do impulsivo no r�dio i
        Pvi = 1/n * sum(abs(V2(i,:)).^2);
        
        if Pvi == 0
            VI(i,:) = V2(i,:) * 0; % Se n�o houver ru�do, zera o vetor
        else
            
            % Normaliza a energia e multiplica por Sigma_i
            VI(i,:) = V2(i,:)/ sqrt(Pvi)* Sigma_i;
        end
    end
    %%
    
    %%
    % RU�DO IMPULSIVO TIPO 2
    
else
    
    % Desvio padr�o para ru�do tipo 2
    Sigma_i = sqrt( K );
    
    % Determina quais r�dios ser�o afetados
    Afected_CR = binornd(1,Pcr,1,m);
    
    % Inicializa matriz complexa de zeros para ru�do impulsivo
    N = zeros(m,n) + 1i * zeros(m,n);
    
    % Para cada r�dio/sensor
    for j = 1 : m
        
        % Calcula a primeira posi��o do ru�do impulsivo via distribui��o exponencial
        position_imp_noise = round(exprnd(Beta));
        
        % Para cada amostra da linha j
        for i = 1 : n
            if i == position_imp_noise
                
                % Gera ru�do normal com m�dia A e desvio B
                X = A + B*randn(1);
                
                % Calcula amplitude lognormal do ru�do impulsivo
                x = 10 .^(X/20);
                
                % Gera fase theta uniformemente entre 0 e 2*pi
                theta = unifrnd(0,2*pi);
                
                % Calcula componentes em fase e quadratura do ru�do impulsivo
                N(j,i) = x .*cos(theta) + 1i* x.*sin(theta);
                
                % Define pr�xima posi��o do ru�do impulsivo, deslocada por vari�vel exponencial
                position_imp_noise = i + round(exprnd(Beta));
            end
        end
    end;
    
    % Remove ru�do dos r�dios n�o afetados
    V2 = diag(Afected_CR) * N;
    
    % Normaliza e aplica fator Sigma_i para cada r�dio
    for i = 1 : m
        
        % Calcula energia do ru�do impulsivo no r�dio i
        Pvi = 1/n * sum(abs(V2(i,:)).^2);
        
        if Pvi == 0 % Se n�o houver ru�do, zera vetor
            VI(i,:) = V2(i,:) * 0;
        else
            
            % Normaliza a energia e multiplica por Sigma_i
            VI(i,:) = V2(i,:)/ sqrt(Pvi)* Sigma_i;
        end
        
    end
    
    % Aplica a probabilidade Pin para ativar o ru�do impulsivo na matriz VI
    VI = binornd(1,Pin)*VI;
end;
%%