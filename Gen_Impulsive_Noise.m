%%
% SIMULA O EFEITO DE RUÍDO IMPULSIVO NAS AMOSTRAS COLETADAS

% Função: Gen_Impulsive_Noise

% Autores: André Antônio dos Anjos e Luis Miguel Alves Borges
%          Universidade Federal de Uberlândia (UFU), Patos de Minas, MG, Brazil
%          Email: {andre.anjos; luis.alves} @ufu.br

% Descrição:
% Simula o efeito do ruído impulsivo sobre sinais complexos coletados
% em m receptores (ou linhas), cada um com n amostras.
% Permite gerar dois tipos de ruído impulsivo (Op_IN_Tipe = 1 ou 2).

% Entradas:
% m - Número de linhas (ex: sensores ou receptores)
% n - Número de amostras por linha
% Nb - Número de bursts (pulsos) de ruído impulsivo (tipo 1)
% Ns - Largura de cada burst de ruído impulsivo (tipo 1)
% K - Relação entre potência do ruído impulsivo e do ruído gaussiano
% Pin - Probabilidade do ruído impulsivo estar presente
% Pcr - Probabilidade de um receptor ser afetado pelo ruído impulsivo
% A, B - Parâmetros do ruído lognormal (tipo 2)
% Beta - Parâmetro da distribuição exponencial para espaçamento (tipo 2)
% Op_IN_Tipe Tipo do ruído impulsivo:
% 1: bursts localizados (pulsos)
% 2: rajadas com espaçamento exponencial
% SNR - Relação sinal-ruído em dB (usado para calcular sigma_i no tipo 1)

% Saídas:
% VI - Matriz m x n de ruído impulsivo complexo
%%

function [VI] = Gen_Impulsive_Noise(m,n,Nb,Ns,K,Pin,Pcr,A,B,Beta,Op_IN_Tipe,SNR)

%%
% RUÍDO IMPULSIVO TIPO 1

if Op_IN_Tipe == 1
    
    % Calcula o desvio padrão do ruído impulsivo baseado na relação SNR e fator K
    Sigma_i = sqrt( K / (10^(SNR/10)));
    
    % Inicializa matriz Z com zeros para indicar posições dos pulsos
    Z = zeros(m,n);
    
    % Para cada linha (para cada rádio ou sensor)
    for u = 1 : m
        
        % Inicializa vetor de zeros para marcar posições de pulso
        g = zeros(1,n);
        
        % Definição dos intervalos entre os pulsos de ruído impulsivo
        t(2 * Nb) = 0;
        
        % Primeiro intervalo aleatório para posicionar o primeiro pulso
        t(1) = randi([0 (floor(n /Nb))]);
        t(2) = t(1) + Ns - 1; % Fim do primeiro pulso
        
        % Se houver múltiplos pulsos (Nb > 1), define os intervalos subsequentes
        if Nb > 1
            for r = 1 : Nb -1
                t(2*r + 1) = t(2*r) + randi([0 (floor(n /Nb))]) + Ns;
                t(2*r + 2) = t(2*r + 1) + Ns - 1;
            end
        end
        
        % Copia os intervalos para z
        z = t;
        
        % Marca as posições dos pulsos no vetor g (com 1)
        for r = 0 : Nb -1
            for w = z(2*r+1) : z(2*r +2)
                g(mod(w,n)+1)= 1;  % Usa módulo para circular no vetor
            end
        end
        
        % Preenche a linha da matriz Z com os pulsos marcados
        Z(u,:) = g;
    end
    
    % Determina quais rádios serão afetados pelo ruído impulsivo, via variável binomial
    for i = 1 : m
        aux_bino(i) = binornd(1,Pcr);
    end
    
    % Remove ruído dos rádios não afetados
    G = diag(aux_bino) * Z;
    
    N = zeros(m,n);
    
    % Gera ruído Gaussiano complexo branco
    N = randn(m,n) +  1i * randn(m,n);
    
    % Aplica ruído somente nas posições indicadas em G
    V2 = N.*G;
    
    % Normaliza e aplica fator Sigma_i para cada rádio
    for i = 1 : m
        
        % Calcula energia do ruído impulsivo no rádio i
        Pvi = 1/n * sum(abs(V2(i,:)).^2);
        
        if Pvi == 0
            VI(i,:) = V2(i,:) * 0; % Se não houver ruído, zera o vetor
        else
            
            % Normaliza a energia e multiplica por Sigma_i
            VI(i,:) = V2(i,:)/ sqrt(Pvi)* Sigma_i;
        end
    end
    %%
    
    %%
    % RUÍDO IMPULSIVO TIPO 2
    
else
    
    % Desvio padrão para ruído tipo 2
    Sigma_i = sqrt( K );
    
    % Determina quais rádios serão afetados
    Afected_CR = binornd(1,Pcr,1,m);
    
    % Inicializa matriz complexa de zeros para ruído impulsivo
    N = zeros(m,n) + 1i * zeros(m,n);
    
    % Para cada rádio/sensor
    for j = 1 : m
        
        % Calcula a primeira posição do ruído impulsivo via distribuição exponencial
        position_imp_noise = round(exprnd(Beta));
        
        % Para cada amostra da linha j
        for i = 1 : n
            if i == position_imp_noise
                
                % Gera ruído normal com média A e desvio B
                X = A + B*randn(1);
                
                % Calcula amplitude lognormal do ruído impulsivo
                x = 10 .^(X/20);
                
                % Gera fase theta uniformemente entre 0 e 2*pi
                theta = unifrnd(0,2*pi);
                
                % Calcula componentes em fase e quadratura do ruído impulsivo
                N(j,i) = x .*cos(theta) + 1i* x.*sin(theta);
                
                % Define próxima posição do ruído impulsivo, deslocada por variável exponencial
                position_imp_noise = i + round(exprnd(Beta));
            end
        end
    end;
    
    % Remove ruído dos rádios não afetados
    V2 = diag(Afected_CR) * N;
    
    % Normaliza e aplica fator Sigma_i para cada rádio
    for i = 1 : m
        
        % Calcula energia do ruído impulsivo no rádio i
        Pvi = 1/n * sum(abs(V2(i,:)).^2);
        
        if Pvi == 0 % Se não houver ruído, zera vetor
            VI(i,:) = V2(i,:) * 0;
        else
            
            % Normaliza a energia e multiplica por Sigma_i
            VI(i,:) = V2(i,:)/ sqrt(Pvi)* Sigma_i;
        end
        
    end
    
    % Aplica a probabilidade Pin para ativar o ruído impulsivo na matriz VI
    VI = binornd(1,Pin)*VI;
end;
%%