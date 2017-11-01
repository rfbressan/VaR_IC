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
  cbind(c("IBovespa", "S&P500", "S&P TSE", "IPSA", "Merval", "IPC"),
        stringsAsFactors = FALSE)

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
## Por testes de auto.arima os modelos ideais de cada serie variam bastante.
## Na media existem mais compenentes MA que AR, portanto um modelo ARMA(1,2)
## parece um bom compromisso entre parcimonia e ajuste.
## Ja as distribuicoes de zt a normal e t-Student nao apresentam bom fit
## A Johnson, GED, NIG, SkewStudent e a Ghyp sao melhores
## Lembrando, o modelo das perdas é ARMA(1,2) e a volatilidade é eGARCH(1,1)
## L_t=mu_t+e_t      mu_t=mu+ar1*mu_t-1+ma1*e_t-1+ma2*e_t-2+e_t
## e_t=sigma_t*z_t   ln(sigma^2_t)=omega+alpha1*z_t-1+gamma1(|z_t-1|-E[|zt-1|])+beta1*ln(sigma^2_t-1)
###################################################################################
# Lfit2 usa método de rugarch
# LEMBRAR: ts contem os retornos e não as perdas!
# Lfit1 <- garchFit(formula = ~arma(1,1)+garch(1,1), 
#                   data = losses[paste0("/", end),],
#                   include.mean = TRUE, 
#                   algorithm = "lbfgsb+nm")

ruspec <- ugarchspec(mean.model = list(armaOrder = c(1,2)),
                     variance.model = list(model = "eGARCH", garchOrder = c(1,1)),
                     distribution.model = "sstd")
garch.models <- assets.tbl[,1:3]
garch.models <- garch.models %>% 
  mutate(garch = map(ts, ~ugarchfit(ruspec, .x[paste0("/", end),], solver = "hybrid")))

# Mostra um sumario para o Ibovespa
show(garch.models$garch[[1]])
# Terminamos ainda com efeito alavancagem a ser modelado para alguns casos
# signbias mostra significancia do efeito negativo
# Talvez um modelo GJR ou APARCH possa resolver. Nao resolveram, melhor eGARCH

plot(-garch.models$ts[[1]][paste0("/", end),], col = "blue")
lines(fitted(garch.models$garch[[1]]), col = "black")

plot(garch.models$garch[[1]])
# Gerar 6 figuras com estes 4 graficos ACF
for(i in 1:dim(garch.models)[1]) {
  jpeg(filename = paste0("artigo-acf-", garch.models$id_name[i], ".jpeg"),
       width = 640, height = 640)
  op <- par(mfrow=c(2,2))
  plot(garch.models$garch[[i]], which = 4)
  plot(garch.models$garch[[i]], which = 5)
  plot(garch.models$garch[[i]], which = 10)
  plot(garch.models$garch[[i]], which = 11)
  par(op)
  dev.off()
}

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

