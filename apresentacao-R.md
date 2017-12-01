



Medidas condicionais de risco com teoria do valor extremo
========================================================
author: Rafael Felipe Bressan
date: 15-11-2017
autosize: true
bibliography: library.bib

Motivação
========================================================

- De acordo com os princípios do acordo de Basileia III, as instituições financeiras supervisionadas pelos Bancos Centrais devem manter *buffers* de capital contra riscos de mercado, crédito, liquidez, entre outros.
- Para riscos de mercado, as duas formas mais usuais de fazer a quantificação destes são os métodos de Valor em Risco - VaR e o *Expected Shortfall* - ES.
- Uma estimação excessiva da medida de risco gerará um excesso de capital em reserva. Custo para a instituição.
- Uma subestimação deste risco pode levar a IF a uma crise de liquidez e eventualmente a insolvência.

Valor em Risco
========================================================

- VaR é um quantil $\alpha$ da distribuição de perdas de um ativo ou portfólio em um determinado período de tempo.
- O método VaR para cálculo de risco de mercado ao qual um portfólio está sujeito foi primeiramente introduzido pelo banco J. P. Morgan em 1995.
- Método original assumia distribuição normal das perdas, correlação constante entre ativos e era calculada de forma incondicional.

Valor em Risco
========================================================

![VaR é um quantil da distribuição de perdas.](apresentacao-R-figure/var-1.png)

Expected Shortfall
========================================================

- ES é o valor esperado das perdas que forem iguais ou maiores que o VaR.
- É calculado como uma média condicional.
- Medida coerente de risco.  @Acerbi2001.

$$
\begin{align*}
VaR_\alpha^t=&\inf\{F_{L_{t+1}} | \mathcal{G}_t(\mathcal{L}) \geq \alpha\}, \\
ES_\alpha^t=&E[L_{t+1} | L_{t+1} > VaR_\alpha^t]
\end{align*}
$$

Expected Shortfall
========================================================

![ES é o valor esperado das perdas, dado que a perda exceda o VaR.](apresentacao-R-figure/ES-1.png)

Fatos estilizados
========================================================

Séries temporais de retornos possuem as seguintes características:

- Ausência de autocorrelação. Componente AR é fraco.

- Grande autocorrelação nos retornos absolutos ou retornos ao quadrado.

- Agrupamento (*Clusters*) de volatilidade.

- Persistência nas autocorrelações dos retornos ao quadrado.

- Distribuição incondicional com caudas longas (leptocúrticas).

- Distribuição condicional com algum grau de leptocurtose.

- Assimetria entre ganhos e perdas.

Fatos estilizados
========================================================

![Retornos do S&P500 de 29/01/2002 até 30/01/2009.](apresentacao-R-figure/sp500ret-1.png)

Modelos GARCH
========================================================

- Modelos GARCH, @Bollerslev1986 lidam com a heteroscedasticidade condicional encontrada nas séries financeiras.

- Propriedades desejáveis: leptocurtose e autocorrelação na variância

- O modelo GARCH exponencial ou eGARCH de @Nelson1991 lida também com o efeito alavancagem.

$$
\begin{align*}
	L_t=&\mu+ \sum_{i=1}^r\phi_i L_{t-i}+ \sum\limits_{j=1}^{s}\theta_j\epsilon_{t-j} +\epsilon_t \\
	\ln(\sigma_t^2)=&\omega+ \sum\limits_{i=1}^{p}(\alpha_i Z_{t-i}+ \gamma_i(|Z_{t-i}|-E|Z_{t-i}|))+ \sum\limits_{j=1}^{q}\beta_j \ln(\sigma_{t-j}^2) \\
	\epsilon_t=&\sigma_t Z_t
\end{align*}
$$

Teoria do valor extremo
========================================================

- O método peaks-over-treshold modela a distribuição dos excessos acima de um determinado limiar.

**Definição** [Distribuição dos excessos]:
	Seja *X* uma variável aleatória com função de distribuição $F \in MDA(H_\xi)$. A distribuição dos excessos sobre um limiar *u* tem a função de Distribuição Generalizada de Pareto - GPD:

$$
\begin{equation*}
  G_{\xi,\beta(u)}(X) = 
  \begin{cases}
    1- \left(1+ \frac{\xi x}{\beta(u)} \right)^{-\frac{1}{\xi}}, & \xi \neq 0,\\
    1-exp\left(-\frac{x}{\beta(u)}\right), & \xi = 0,\\
  \end{cases}
\end{equation*}
$$

Os parâmetros $\xi$ e $\beta$ são conhecidos respectivamente como parâmetros de forma e escala da distribuição. 

Teoria do valor extremo
========================================================

![Método POT para estimar caudas de distribuições.](apresentacao-R-figure/pot-1.png)

Modelo eGARCH-EVT
========================================================

- Seguiremos os passos propostos por @McNeil2000.
- Nosso modelo completo para as medidas de risco $VaR_\alpha$ e $ES_\alpha$ condicionais dada a distribuição de perdas $L_t$ de um ativo será, portanto:

$$
\begin{align*}
L_t=&\mu+ \phi_1 L_{t-1}+ \epsilon_t \\
\epsilon_t=&\sigma_t Z_t\\
\ln(\sigma_t^2)=&\omega+ \sum_{i=1}^{2}(\alpha_i Z_{t-i}+ \gamma_i(|Z_{t-i}|-E|Z_{t-i}|))+ \beta_1 \ln(\sigma_{t-1}^2) \\
Z_t\sim &\mathcal{D}(0,1) \text{ e } \mathcal{D} \in MDA(H_\xi)
\end{align*}
$$

VaR e ES parametrizados
========================================================

- Nossas medidas de risco neste modelo condicional serão:

$$
\begin{align*}
VaR_\alpha^t=&\mu_{t+1}+\sigma_{t+1}z_\alpha, \\
ES_\alpha^t=&\mu_{t+1}+\sigma_{t+1}E[Z | Z>z_\alpha]
\end{align*}
$$

onde $z_\alpha$ é o quantil $\alpha$ das inovações *Z*.

VaR condicional estimado
========================================================

\begincols

  \begincol{.65\textwidth}
  
  ![VaR condicional para o S\&P500.](artigo-apresentacao-var.jpeg)
  
  \endcol

  \begincol{.31\textwidth}
  
  \input{artigo-apresentacao-tabela.tex}
  
  \endcol
  
\endcols

Referências