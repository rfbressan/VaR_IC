## Artigo a ser proposto na disciplina MEFCA. Baseado no método de McNeil2000
## Filtar os retornos diários de vários índices primeiramente por um modelo 
## ARMA-GARCH para em seguida ajustar uma distribuição GPD aos resíduos.
## Com esta modelagem é possível estimar os valores de VaR e ES
# Indices de bolsas utilizados
# BVSP - Bovespa Brasil
# MERV - Merval Argentina
# IPSA - IPSA Chile
# MXX - IPC México
# GSPC - SP500 EUA
# GSPTSE - SP/TSX Canada
# Dados entre 31/08/2005 a 31/08/2017

# Inicio ------------------------------------------------------------------

#library(fGarch)
library(fExtremes)
library(fBasics)
library(QRM)
library(rugarch)
library(timeSeries)
library(xts)
library(PerformanceAnalytics)
library(xtable)
library(tidyverse)
library(broom)
library(purrr)
library(gridExtra)
library(ggplot2)
library(CADFtest) # Teste Dickey-Fuller

start <- as.Date("2005-08-31")
end <- as.Date("2016-08-31")
backstart <- as.Date("2016-09-01")

list.returns <- function(asset) {
  tb <- read.csv(paste0("artigo-", asset, ".csv"), stringsAsFactors = FALSE)
  prices <- as.xts(read.zoo(tb, format = "%Y-%m-%d", FUN = as.Date))
  return(100*na.omit(Return.calculate(prices$Adj.Close, method = "log")))
}
# Gera um tible com uma coluna com o codigo do ativo, a serie de retornos - ts e
# o nome do indice - id_name
assets <- c("BVSP", "GSPC", "GSPTSE", "IPSA", "MERV", "MXX")
lista <- lapply(assets, list.returns)
names(lista) <- assets
assets.tbl <- enframe(lista) %>% 
  cbind(c("IBovespa", "S&P500", "S&P TSE", "IPSA", "Merval", "IPC"),
        stringsAsFactors = FALSE)

colnames(assets.tbl) <- c("indice", "ts", "id_name")
# Qual o menor numero de dias entre inicio do periodo de backtesting ate o final da serie?
outsample = enframe(map_dbl(lista, ~ndays(.x[paste0(backstart, "/")])))
names(outsample) <- c("indice", "out")
# Remove a variavel lista que agora e desnecessaria
rm(lista)

# Estatisticas descritivas retornos-----------------------------------------
assets.tbl <- assets.tbl %>% 
  mutate(media = map_dbl(ts, ~ mean(.x)),
         mediana = map_dbl(ts, ~median(.x)),
         maximo = map_dbl(ts, ~max(.x)),
         minimo = map_dbl(ts, ~min(.x)),
         desvp = map_dbl(ts, ~sd(.x)),
         assim = map_dbl(ts, ~skewness(.x)),
         curtose = map_dbl(ts, ~kurtosis(.x)),
         jarquebera = map(ts, ~jarqueberaTest(as.timeSeries(.x))),
         ljungbox = map(ts, ~Box.test(.x, lag = 10, type = "Ljung-Box")),
         ljungbox2 = map(ts, ~Box.test(.x^2, lag = 10, type = "Ljung-Box"))) %>% 
  mutate(jbstat = map_dbl(jarquebera, ~as.numeric(.x@test$statistic)),
         jbpvalue = map_dbl(jarquebera, ~as.numeric(.x@test$p.value)),
         lbstat = map_dbl(ljungbox, ~as.numeric(.x$statistic)),
         lbpvalue = map_dbl(ljungbox, ~as.numeric(.x$p.value)),
         lb2stat = map_dbl(ljungbox2, ~as.numeric(.x$statistic)),
         lb2pvalue = map_dbl(ljungbox2, ~as.numeric(.x$p.value)),
         nobs = map_int(ts, ~as.integer(length(.x))))

# Transpor esta tabela para apresentacao
df.descritivas <- data.frame(t(assets.tbl[,-c(1, 2, 3, 11, 12, 13)]))
colnames(df.descritivas) <- assets.tbl$id_name
row.names(df.descritivas) <- c("Média", "Mediana", "Máximo", "Mínimo", "Desvp", "Assimetria", "Curtose",
                               "Jarque-Bera", "p-valor", "Q(10)", "p-valo", "Q^2(10)", "p-val",
                               "N.obs")
# Cria o xtable
tab1 <- xtable(df.descritivas, caption = "Estatísticas descritivas dos retornos 
               (amostra completa de 31/08/2005 a 30/08/2017).",
               digits = 5,
               label = "tab:descritivas",
               auto = TRUE)
print.xtable(tab1, 
             file = "artigo-tab-descritivas.tex",
             caption.placement = "top",
             table.placement = "H")
# Graficos retornos------------------------------------------------------
list.plot <- lapply(seq_along(assets.tbl$indice), 
                    function(x) {autoplot(assets.tbl$ts[[x]])+
                        labs(x = "", y = "", title = paste(assets.tbl$id_name[[x]], "retornos"))}) 

grid.arrange(grobs = list.plot)

## Teste de estacionariedade das series de retornos
# Teste Dickey-Fuller encontrado no pacote CADFtest
# Delta_y = mu + theta*t + delta*y_t-1
adf <- lapply(assets.tbl$ts, CADFtest)
names(adf) <- assets.tbl$indice
adf <- enframe(adf)
adf <- adf %>%
  mutate(teste = map(value, ~tidy(.x)))
adf$value <- NULL
adf <- unnest(adf)
adf

# Teste para graficos QQ normal
op <- par(mfrow = c(3,2))
for(i in 1:dim(assets.tbl)[1]){
  qqnormPlot(as.timeSeries(assets.tbl$ts[[i]]), labels = FALSE, title = FALSE, mtext = FALSE,
             main = assets.tbl$id_name[[i]],
             xlab = "Normal",
             ylab = "Amostra")
}
par(op)


# Modelo eGARCH -----------------------------------------------------------
## Ja as distribuicoes de zt a normal e t-Student nao apresentam bom fit
## A Johnson, GED, NIG, SkewStudent e a Ghyp sao melhores
## Lembrando, o modelo das perdas é AR(1) e a volatilidade é eGARCH(2,1)
## L_t=mu_t+e_t      mu_t=mu+ar1*mu_t-1+e_t
## e_t=sigma_t*z_t   ln(sigma^2_t)=omega+alpha1*z_t-1+gamma1(|z_t-1|-E[|zt-1|])+beta1*ln(sigma^2_t-1)
## LEMBRAR: ts contem os retornos e não as perdas!

ruspec <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "eGARCH", garchOrder = c(2,1)),
                     distribution.model = "sstd")

## Modelando as PERDAS!! parametro eh loss para a funcao ugarchfit
# Deixamos um numero outsample para fazer o backtesting. Diferente para cada ativo
garch.models <- assets.tbl[,1:3] %>% 
  inner_join(outsample, by = "indice") %>% 
  mutate(loss = map(ts, ~-.x),
         garch = map2(loss, out, ~ugarchfit(ruspec, 
                                            .x,
                                            out.sample = .y,
                                            solver = "hybrid")),
         ts = NULL)

# Mostra a convergencia para cada modelo
lapply(garch.models$garch, convergence)

# Sumarios dos modelos GARCH
show(garch.models$garch[[1]])

# Terminamos ainda com efeito alavancagem a ser modelado para alguns casos
# signbias mostra significancia do efeito negativo
# Talvez um modelo GJR ou APARCH possa resolver. Nao resolveram, melhor eGARCH(2,1)

# Gerar 6 figuras com estes 4 graficos ACF
for(i in 1:dim(garch.models)[1]) {
  jpeg(filename = paste0("artigo-acf-", garch.models$id_name[i], ".jpeg"),
       width = 800, height = 800, quality = 100)
  op <- par(mfrow=c(2,2))
  plot(garch.models$garch[[i]], which = 4)
  plot(garch.models$garch[[i]], which = 5)
  plot(garch.models$garch[[i]], which = 10)
  plot(garch.models$garch[[i]], which = 11)
  par(op)
  dev.off()
}
file.rename(c("artigo-acf-S&P500.jpeg", "artigo-acf-S&P TSE.jpeg"), 
            c("artigo-acf-SP500.jpeg", "artigo-acf-SP-TSE.jpeg"))
# Cria tabela com os parametros estimados e seus p-valores
vec_coef <- function(model) {
  vec(t(model@fit$robust.matcoef[, c(1,4)]))
}

mat <- do.call(cbind,
                 lapply(garch.models$garch, vec_coef))
param <- c("$\\mu$", "",
           "$\\phi_1$", "",
           "$\\omega$", "",
           "$\\alpha_1$", "",
           "$\\alpha_2$", "",
           "$\\beta_1$", "",
           "$\\gamma_1$", "",
           "$\\gamma_2$", "",
           "$\\zeta$", "",
           "$\\nu$", "")
garchcoef <- data.frame(par = param, mat)

colnames(garchcoef) <- c("Parâmetros", garch.models$id_name)
# Xtable
tab2 <- xtable(garchcoef, caption = "Par\\^ametros estimados do modelo eGARCH. Valores p apresentados 
               de acordo com erros padrão robustos. (amostra de trabalho entre 31/08/2005 a 31/08/2016).",
               digits = 5,
               label = "tab:garchcoef",
               auto = TRUE)
print.xtable(tab2, 
             file = "artigo-tab-garchcoef.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

# Modelo EVT para os residuos padronizados --------------------------------

## Os residuos podem ser retirados atraves do metodo residuals com a opcao "standardize=T"
## a coluna rq ira guardar os valores de zq (quantile) e sq (shortfall)
evt.models <- garch.models %>% 
  mutate(resid_z = map(garch, ~coredata(residuals(.x, standardize = TRUE))),
         Nu = map(resid_z, ~sum(.x > quantile(.x, 0.95))),
         gpdfit = map(resid_z, ~gpdFit(.x, u = quantile(.x, 0.95))),
         u = map_dbl(gpdfit, ~.x@parameter$u),
         xi = map_dbl(gpdfit, ~.x@fit$par.ests[1]),
         beta = map_dbl(gpdfit, ~.x@fit$par.ests[2]),
         xi_se = map_dbl(gpdfit, ~.x@fit$par.ses[1]),
         beta_se = map_dbl(gpdfit, ~.x@fit$par.ses[2]),
         rq = map(gpdfit, ~gpdRiskMeasures(.x, prob = c(0.975, 0.99))),
         loss = NULL,
         garch = NULL)
## Teste com o pacote evir
testez <- coredata(residuals(garch.models$garch[[1]], standardize = TRUE))
teste_evir <- gpd(testez, threshold = quantile(testez, 0.95)) # Mesmos valores do fExtremes
meplot(-testez)
shape(testez, models = 10, start = 90, end = 150)
## Teste com o pacote evd
teste_evd <- fpot(testez, quantile(testez, 0.95))
fitted.values(teste_evd)
std.errors(teste_evd)  # Iguais ao fExtremes e evir
mrlplot(-testez)

## Gráficos para analisar a qualidade do gpdFit
op <- par(mfrow=c(2,2))
plot(evt.models$gpdfit[[1]], which='all')
par(op)

# Reconstruindo o VaR e o ES condicionais ---------------------------------

## VaR: xq_t = mu_t+1 + sigma_t+1*zq
## ES: Sq_t = mu_t+1 + sigma_t+1*sq

# Teste para in sample VaR
zq <- evt.models$rq[[2]]$quantile[1]
xq <- fitted(garch.models$garch[[2]])+sigma(garch.models$garch[[2]])*zq
xql <- xq[-length(xq)]
linsample <- garch.models$loss[[2]][paste0("/", end)][-1]
plot(linsample)
lines(xql, col = "red")
varex <- sum(coredata(linsample) > coredata(xql))
varex/length(xql)*100

# E por fim calcula as medidas de risco para os residuos zt
risks <- gpdRiskMeasures(evtfit, prob = 0.99) # Medidas sem intervalo de conf.
tail <- gpdTailPlot(evtfit)
varci <- gpdQPlot(tail)
esci <- gpdSfallPlot(tail)
risktable <- rbind(varci, esci)
dimnames(risktable) <- list(c("$z_{.99}$", "$s_{.99}$"), c("Inf", "Estimativa", "Sup"))
print.xtable(xtable(risktable, caption = "Valores de $z_{.99}$ e $s_{.99}$ encontrados e seus
                    respectivos intervalos de confian\\c ca a 95\\%",
                    label = "tab:tabevtAAPL2", align = c("r", "r", "c", "r")),
             file = "..\\tables\\tabevtAAPL2.tex", sanitize.text.function = function(x) {x})

