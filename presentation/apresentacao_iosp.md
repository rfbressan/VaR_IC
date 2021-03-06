---
title: "Aplicando a teoria do valor extremo no cálculo de risco de índices setoriais da Bovespa"
author: "Rafael Felipe Bressan"
date: "2018-08-17"
output:
  iosp::ioslides_plus:
    footer: "Grupo de Pesquisa em Economia Aplicada - GEA"
    widescreen: false
    logo: ../artigo/figs/Marca_Udesc.jpeg
    css: ['css/box.css', 'css/lecture.css']
    smaller: true
    keep_md: true
---



## Agenda

> - Motivação
> - Valor em Risco - VaR
> - Teoria do Valor Extremo - EVT
> - Modelo EGARCH
> - Índices Utilizados
> - Resultados

# Introdução

## Motivação 

> - Princípios de Basileia 
> - Instituições Financeiras devem manter reservas de capital contra riscos de mercado, crédito, liquidez, entre outros
> - Riscos de mercado,  Valor em Risco - VaR e *Expected Shortfall* - ES.
> - Estimação excessiva do risco gera excesso de capital em reserva. Custo para a instituição.
> - Subestimação deste risco pode levar a IF a uma crise de liquidez e eventualmente a insolvência.

## Valor em Risco

- VaR é um quantil $\alpha$ da distribuição de perdas de um ativo ou portfólio em um determinado período de tempo.


<div id = "var" class = " ">
<img src="apresentacao_iosp_files/figure-html/var-1.png" style="display: block; margin: auto;" />
</div>

# Fundamentação Teórica

## Teoria do valor extremo

- O método peaks-over-treshold modela a distribuição dos excessos acima de um determinado limiar.

**Definição** [Distribuição dos excessos]:
	Seja *X* uma variável aleatória com função de distribuição $F \in MDA(H_\xi)$. A distribuição dos excessos sobre um limiar *u* tem a função de Distribuição Generalizada de Pareto - GPD:

\begin{equation*}
  G_{\xi,\beta(u)}(X) = 
  \begin{cases}
    1- \left(1+ \frac{\xi x}{\beta(u)} \right)^{-\frac{1}{\xi}}, & \xi \neq 0,\\
    1-exp\left(-\frac{x}{\beta(u)}\right), & \xi = 0,\\
  \end{cases}
\end{equation*}
Os parâmetros $\xi$ e $\beta$ são conhecidos respectivamente como parâmetros de forma e escala da distribuição. 

## Filtro EGARCH

- O modelo GARCH exponencial ou EGARCH de @Nelson1991 lida também com o efeito alavancagem, além da heterocedasticidade

\begin{align*}
	L_t=&\mu+ \sum_{i=1}^r\phi_i L_{t-i}+ \sum\limits_{j=1}^{s}\theta_j\epsilon_{t-j} +\epsilon_t \\
	\ln(\sigma_t^2)=&\omega+ \sum\limits_{i=1}^{p}(\alpha_i Z_{t-i}+ \gamma_i(|Z_{t-i}|-E|Z_{t-i}|))+ \sum\limits_{j=1}^{q}\beta_j \ln(\sigma_{t-j}^2) \\
	\epsilon_t=&\sigma_t Z_t
\end{align*}

# Modelo

## Modelo eGARCH-EVT

- Seguiremos os passos propostos por @McNeil2000.
- Modelo completo para as medidas de risco $VaR_\alpha$ e $ES_\alpha$ condicionais dada a distribuição de perdas $L_t$ de um ativo será, portanto:

\begin{align*}
L_t=&\mu+ \phi_1 L_{t-1}+ \epsilon_t \\
\epsilon_t=&\sigma_t Z_t\\
\ln(\sigma_t^2)=&\omega+ \sum_{i=1}^{2}(\alpha_i Z_{t-i}+ \gamma_i(|Z_{t-i}|-E|Z_{t-i}|))+ \beta_1 \ln(\sigma_{t-1}^2) \\
Z_t\sim &\mathcal{D}(0,1) \text{ e } \mathcal{D} \in MDA(H_\xi)
\end{align*}

## VaR e ES parametrizados

- Medidas de risco neste modelo condicional serão:

\begin{align*}
VaR_\alpha^t=&\mu_{t+1}+\sigma_{t+1}z_\alpha, \\
ES_\alpha^t=&\mu_{t+1}+\sigma_{t+1}E[Z | Z>z_\alpha]
\end{align*}
onde $z_\alpha$ é o quantil $\alpha$ das inovações *Z*.

## Comparação EVT condicioanl _versus_ Normal incondicional

<!-- ![S&P500](./figs/artigo-sp500backtest) -->

# Dados utilizados

## Seis índices de ações

<!-- ![Retornos logarítimicos dos índices de ações](./figs/artigo-retornos_2_3.pdf) -->

## Estatísticas descritivas dos retornos

## Descritivas pós filtragem

<!-- \include(tables/artigo-tab-garchstats) -->

# Resultados

## Teste fora da amostra

- Período de _backtest_ vai de 01/01/2009 a 30/08/2017
- Avaliação dos modelos através de três testes estatístios
- Cobertura incondicional de Kupiec
- Independência entre violações de Christoffersen & Pelletier
- Conjunto de confiança de modelos de Hansen

## Teste de Kupiec

<!-- \include(tables/artigo-tab-varviol) -->

## Sumário de rejeições de modelo bem especificado

<!-- \include(tables/artigo-tab-vartest_suma) -->

## MCS de Hansen

<!-- \include(tables/artigo-tab-mcs) -->

## Referências
