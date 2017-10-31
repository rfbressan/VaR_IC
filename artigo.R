###############################################################################
## Artigo a ser proposto na disciplina MEFCA. Baseado no método de McNeil2000
## Filtar os retornos diários de vários índices primeiramente por um modelo 
## ARMA-GARCH para em seguida ajustar uma distribuição GPD aos resíduos.
## Com esta modelagem é possível estimar os valores de VaR e ES
###############################################################################
# Indices de bolsas utilizados
# BVSP - Bovespa Brasil
# MERV - Merval Argentina
# IPSA - IPSA Chile
# MXX - IPC México
# GSPC - SP500 EUA
# GSPTSE - SP/TSX Canada
# Dados entre 31/08/2005 a 31/08/2017

library(fGarch)
library(fExtremes)
library(QRM)
library(rugarch)
library(timeSeries)
library(xts)
library(PerformanceAnalytics)
library(xtable)
library(tidyverse)
library(purrr)
library(gridExtra)
library(ggplot2)

start <- as.Date("2005-08-31")
end <- as.Date("2016-08-31")
backstart <- as.Date("2016-09-01")

list.returns <- function(asset) {
  tb <- read.csv(paste0("artigo-", asset, ".csv"), stringsAsFactors = FALSE)
  #tb <- tb[-which(tb$Adj.Close == "null"),] # remove linhas com "null"
  #tb[,2:7] <- lapply(tb[,2:7], as.numeric)
  prices <- as.xts(read.zoo(tb, format = "%Y-%m-%d", FUN = as.Date))
  return(100*na.omit(Return.calculate(prices$Adj.Close, method = "log")))
  
}
# Gera um tible com uma coluna com o codigo do ativo, a serie de retornos - ts e
# o nome do indice - id_name
assets <- c("BVSP", "GSPC", "GSPTSE", "IPSA", "MERV", "MXX")
lista <- lapply(assets, list.returns)
names(lista) <- assets
assets.tbl <- enframe(lista) %>% 
  cbind(c("IBovespa", "S&P500", "S&P TSE", "IPSA", "Merval", "IPC"))

colnames(assets.tbl) <- c("indice", "ts", "id_name")
# Remove a variavel lista que agora e desnecessaria
rm(lista)

# Estatisticas descritivas ------------------------------------------------
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
               (amostra completa de 31/08/2005 a 30/08/2017.",
               digits = 5,
               label = "tab:descritivas",
               auto = TRUE)
print.xtable(tab1, 
             file = "artigo-tab-descritivas.tex",
             caption.placement = "top",
             table.placement = "ht")
# Graficos ----------------------------------------------------------------
list.plot <- lapply(seq_along(assets.tbl$indice), 
                    function(x) {autoplot(assets.tbl$ts[[x]])+
                        labs(x = "", y = "", title = paste(assets.tbl$id_name[[x]], "retornos"))}) 

grid.arrange(grobs = list.plot)

# Teste para graficos QQ normal
op <- par(mfrow = c(3,2))
for(i in 1:dim(assets.tbl)[1]){
  qqnormPlot(as.timeSeries(assets.tbl$ts[[i]]), labels = FALSE, title = FALSE, mtext = FALSE,
             main = assets.tbl$id_name[[i]],
             xlab = "Normal",
             ylab = "Amostra")
}
par(op)
###################################################################################
## Lembrando, o modelo das perdas é ARMA(1,1) e a volatilidade é EGARCH(1,1)
## L_t=mu_t+e_t      mu_t=mu+phi*mu_t-1+theta*e_t-1
## e_t=sigma_t*z_t   sigma^2_t=alpha0+alpha1*e^2_t-1+beta*sigma^2t-1
## O modelo retorno 5 parametros:
## mu = valor do intercepto da equacao da media, nosso modelo=0 mas pode não ser
## ar1 = phi do modelo
## ma1 = theta do modelo
## omega = alpha0 do modelo
## alpha1 = alpha1 do modelo
## beta1 = beta do modelo
###################################################################################

# Lfit1 usa método de fGarch
# Lfit2 usa método de rugarch
# LEMBRAR: ts contem os retornos e não as perdas!
Lfit1 <- garchFit(formula = ~arma(1,1)+garch(1,1), 
                  data = losses[paste0("/", end),],
                  include.mean = TRUE, 
                  algorithm = "lbfgsb+nm")

ruspec <- ugarchspec(mean.model = list(armaOrder = c(1,1)),
                     variance.model = list(model = "eGARCH", garchOrder = c(1,1)))
Lfit2 <- ugarchfit(ruspec, losses[paste0("/", end),], solver = "hybrid")
show(Lfit2)
op <- par(mfrow=c(2,3))
plot(Lfit2, which = 1)
plot(Lfit2, which = 4)
plot(Lfit2, which = 5)
plot(Lfit2, which = 9)
plot(Lfit2, which = 10)
plot(Lfit2, which = 11)
par(op)

matcoef <- Lfit2@fit$matcoef
dimnames(matcoef) <- list(c("$\\mu$", "$\\phi_1$", "$\\theta_1$", 
                            "$\\omega$", "$\\alpha_1$", "$\\beta_1$", "$\\gamma_1$"),
                          c("Estimativa", "Erro Padr\\~ao", "Valor t", "Pr(>|t|)"))

print.xtable(xtable(matcoef, caption = "Par\\^ametros estimados para o modelo ARMA-GARCH.",
                    label = "tab:artigoarma", digits = 4),
             file = "tabartigoarma.tex", sanitize.text.function = function(x) {x})

#########################################################################################
## Modelo EVT para os residuos padronizados
## Os residuos podem ser retirados atraves do metodo residuals com a opcao "standardize=T"
##

zt <- as.timeSeries(residuals(Lfit2, standardize = TRUE))

mrlPlot(zt) # Mean Residual Life para escolher treshold u
# u por volta de 1.5% parece ser um valor adequado
# como o quantil 95% eh 1.62% vamos manter o quantil que eh o padrao
# da funcao gpdFit

# Contagem do numero de excessos apenas para verificar se eh um numero
# razoavelmente alto
sum(zt>quantile(zt, 0.95)) # 97, OK

evtfit <- gpdFit(zt, u = quantile(zt, 0.95))
op <- par(mfrow=c(2,2))
plot(evtfit, which='all')
par(op)

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

