# Sensoriamento_Espectral_PFC

Script para simulação e comparação de diferentes técnicas de sensoriamento espectral

Autores: André Antônio dos Anjos e Luis Miguel Alves Borges

Universidade Federal de Uberlândia (UFU), Patos de Minas, MG, Brasil

Email: {andre.anjos; luis.alves} @ufu.br

Este conjunto de scripts MATLAB realiza a simulação do desempenho de técnicas de sensoriamento espectral sob diferentes condições de canal e ruído, gerando curvas ROC (simuladas e teóricas) para avaliar a probabilidade de detecção (Pd) e de falso alarme (Pfa), além de calcular a área sob a curva ROC (AUC) como métrica de desempenho.

O script principal, denominado Codigo_PFC_Final, é configurado em três seções:

Configurações de parâmetros: define o número de experimentos de Monte Carlo, o número total de amostras, SNR média em dB, fator de amostragem, número de símbolos M-PSK, bins para estimar a PDF de fase, tipo de canal, técnica de sensoriamento espectral, incerteza do ruído, erro na frequência central e o número de pontos do limiar para gerar a curva ROC.

Parâmetros default dos canais: permite configurar detalhes específicos de cada modelo de canal, como Rayleigh, Rice, Nakagami-m, Alpha-Mu, Kappa-Mu e Eta-Mu.

Parâmetros do ruído impulsivo: possibilita uma modelagem detalhada deste tipo de ruído, incluindo número de RCs simulados, fator multiplicativo da variância térmica (K), probabilidades de ocorrência, amplitude média, variância dos pulsos e intervalo médio. Para desativar o ruído impulsivo, basta configura-se INParam.K = 0 e INParam.Pin = 0.

Após definir os parâmetros desejados, execute o script Codigo_PFC_Final no MATLAB. O programa gerará automaticamente:

Arquivos .dat com valores da AUC.

Arquivos .fig contendo os gráficos das curvas ROC.

Arquivos .dat adicionais com vetores de Pd, Pfa, thresholds e curvas teóricas.

Um .mat consolidando todos os parâmetros usados, garantindo reprodutibilidade.

Os nomes dos arquivos são gerados dinamicamente, incorporando os principais parâmetros (ex.: N10000_SNR_db_-10_MCe_20000_AUC.dat), o que facilita a organização dos experimentos.

Para executar diferentes cenários de simulação, o usuário pode alterar o tipo de canal por meio do parâmetro Channel, modificar a SNR em dB, escolher diferentes técnicas de sensoriamento espectral ajustando Sepec_Sens_Technique, ou ainda alterar N_pts para gerar curvas ROC com maior resolução.
