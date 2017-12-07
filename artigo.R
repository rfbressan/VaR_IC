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
list.of.packages <- c("fExtremes", "rugarch", "xts", "PerformanceAnalytics", "xtable", "tidyverse", 
                      "broom", "purrr", "gridExtra", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#library(fGarch)
library(fExtremes)
#library(fBasics)
#library(QRM)
library(rugarch)
#library(timeSeries)
library(xts)
library(PerformanceAnalytics)
library(xtable)
library(tidyverse)
library(broom)
library(purrr)
library(gridExtra)
library(ggplot2)
#library(CADFtest) # Teste Dickey-Fuller
source("artigo_fun.R") # Carrega a funcao roll_fit para fazer o backtest


start <- as.Date("2005-08-31")
end <- as.Date("2013-08-31")
backstart <- as.Date("2013-09-01")

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
outsample <-  enframe(map_dbl(lista, ~ndays(.x[paste0(backstart, "/")])))
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
         garch_fit = map2(loss, out, ~ugarchfit(ruspec, 
                                            .x,
                                            out.sample = .y,
                                            solver = "hybrid")),
         ts = NULL)

# Mostra a convergencia para cada modelo
lapply(garch.models$garch_fit, convergence)

# Sumarios dos modelos GARCH
show(garch.models$garch_fit[[1]])

## Para obter os valores de mu_t, sigma_t e z_t fora da amostra
## se utiliza o metodo ugarchfilter com os parametros fixados nos 
## valores estimados pelo ugarchfit
garch.filtered <- garch.models %>% 
  mutate(n_old = map_int(garch_fit, ~.x@model$modeldata$T),
         f_params = map(garch_fit, ~as.list(coef(.x))),
         out_spec = map(f_params, ~ugarchspec(mean.model = list(armaOrder = c(1,0)),
                                              variance.model = list(model = "eGARCH", garchOrder = c(2,1)),
                                              distribution.model = "sstd",
                                              fixed.pars = .x))) %>% 
  mutate(garch_filter = pmap(., ~ugarchfilter(..8,
                                              ..4,
                                              n.old = ..6)))

# Teste para verificar os residuos padronizados sao os mesmos inicialmente
head(rugarch::residuals(garch.filtered$garch_fit[[2]], standardize = TRUE))
head(rugarch::residuals(garch.filtered$garch_filter[[2]], standardize = TRUE))
# Agora garch.filtered contem tudo sobre os modelos Garch, nao precisamos mais de garch.models
rm(garch.models)

# Terminamos ainda com efeito alavancagem a ser modelado para alguns casos
# signbias mostra significancia do efeito negativo
# Talvez um modelo GJR ou APARCH possa resolver. Nao resolveram, melhor eGARCH(2,1)

# Gerar 6 figuras com estes 4 graficos ACF
for(i in 1:dim(garch.filtered)[1]) {
  jpeg(filename = paste0("artigo-acf-", garch.filtered$id_name[i], ".jpeg"),
       width = 800, height = 800, quality = 100)
  op <- par(mfrow=c(2,2))
  plot(garch.filtered$garch_fit[[i]], which = 4)
  plot(garch.filtered$garch_fit[[i]], which = 5)
  plot(garch.filtered$garch_fit[[i]], which = 10)
  plot(garch.filtered$garch_fit[[i]], which = 11)
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
                 lapply(garch.filtered$garch_fit, vec_coef))
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

colnames(garchcoef) <- c("Parâmetros", garch.filtered$id_name)
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
## os valores zq (quantile) e sq (shortfall)
evt.models <- garch.filtered %>% 
  transmute(indice = indice,
            id_name = id_name,
            loss = loss,
            out = out,
            n_old = n_old,
            resid_z = map(garch_filter, ~coredata(residuals(.x, standardize = TRUE))),
            mut = map(garch_filter, ~coredata(fitted(.x))),
            sigmat = map(garch_filter, ~coredata(sigma(.x))),
            Nu = map_int(resid_z, ~sum(.x > quantile(.x, 0.95))),
            gpdfit = map(resid_z, ~gpdFit(.x, u = quantile(.x, 0.95))),
            u = map_dbl(gpdfit, ~.x@parameter$u),
            xi = map_dbl(gpdfit, ~.x@fit$par.ests[1]),
            beta = map_dbl(gpdfit, ~.x@fit$par.ests[2]),
            xi_se = map_dbl(gpdfit, ~.x@fit$par.ses[1]),
            beta_se = map_dbl(gpdfit, ~.x@fit$par.ses[2]),
            zq975= map_dbl(gpdfit, ~gpdRiskMeasures(.x, prob = 0.975)$quantile),
            zq990= map_dbl(gpdfit, ~gpdRiskMeasures(.x, prob = 0.990)$quantile),
            sq975= map_dbl(gpdfit, ~gpdRiskMeasures(.x, prob = 0.975)$shortfall),
            sq990= map_dbl(gpdfit, ~gpdRiskMeasures(.x, prob = 0.990)$shortfall))
# ## Teste com o pacote evir
# testez <- coredata(residuals(garch.models$garch_fit[[1]], standardize = TRUE))
# teste_evir <- gpd(testez, threshold = quantile(testez, 0.95)) # Mesmos valores do fExtremes
# meplot(-testez)
# shape(testez, models = 10, start = 90, end = 150)
# ## Teste com o pacote evd
# teste_evd <- fpot(testez, quantile(testez, 0.95))
# fitted.values(teste_evd)
# std.errors(teste_evd)  # Iguais ao fExtremes e evir
# mrlplot(-testez)

## Gráficos para analisar a qualidade do gpdFit
op <- par(mfrow=c(2,2))
plot(evt.models$gpdfit[[1]], which='all')
par(op)

# Reconstruindo o VaR e o ES condicionais ---------------------------------

## VaR: xq_t = mu_t+1 + sigma_t+1*zq
## ES: Sq_t = mu_t+1 + sigma_t+1*sq
riskmeasures <- evt.models %>% 
  transmute(indice = indice,
            id_name = id_name,
            loss = loss,
            out = out,
            n_old = n_old,
            VaR975 = pmap(., ~(..7+..8*..16)), ## VaR = mu_t+1 + sigma_t+1*zq
            VaR990 = pmap(., ~(..7+..8*..17)),
            ES975 = pmap(., ~(..7+..8*..18)),  ## ES = mu_t+1 + sigma_t+1*sq
            ES990 = pmap(., ~(..7+..8*..19))) %>% 
  mutate(out_VaR975 = pmap(., ~..6[c((..5+1):(..5+..4-1))]), #out_VaR = VaR[(n_old+1):(n_old+out-1)]
         out_VaR990 = pmap(., ~..7[c((..5+1):(..5+..4-1))]),
         out_ES975 = pmap(., ~..8[c((..5+1):(..5+..4-1))]),  #out_ES = ES[(n_old+1):(n_old+out-1)]
         out_ES990 = pmap(., ~..9[c((..5+1):(..5+..4-1))]),
         out_loss = map2(loss, n_old, ~coredata(.x[-c(1:(.y+1))])[, 1, drop = TRUE])) #out_los = loss[-c(1:n_old+1)]
# Plotando os valores fora da amostra
plot_risks <- function(loss, VaR, ES, id_name) {
  tindex <- 1:length(loss)
  df <- data.frame(x = tindex, loss = loss, VaR = VaR, ES = ES)
  plot <- ggplot(df, aes(x = x))+
    geom_line(aes(y = loss), color = "black")+
    geom_line(aes(y = VaR), color = "red")+
    geom_line(aes(y = ES), color = "darkgreen")+
    labs(x = "", y = "", title = id_name)
  return(plot)
}
VaR_plots <- riskmeasures %>% 
  transmute(VaR975_plot = pmap(., ~plot_risks(..14, ..10, ..12, ..2)),
            VaR990_plot = pmap(., ~plot_risks(..14, ..11, ..13, ..2)))

grid.arrange(grobs = VaR_plots$VaR975_plot)
grid.arrange(grobs = VaR_plots$VaR990_plot)

# Teste APENAS PARA FORA DA AMOSTRA 
## Com os valores de n_old, out e as perdas, calcular a quantidade de violacoes do VaR 
## e colocar em um tibble para cada ativo
#VaR para uma normal incondicional
out_loss <- riskmeasures$out_loss[[1]]
varnorm <- qnorm(0.990,
                 mean(coredata(garch.filtered$loss[[1]][paste0("/", end)])),
                 sd(coredata(garch.filtered$loss[[1]][paste0("/", end)])))
varex <- sum(out_loss > riskmeasures$out_VaR990[[1]])
varex/length(out_loss)*100
varexnorm <- sum(out_loss > varnorm)
varexnorm/length(out_loss)*100
dfex <- data.frame(modelo = c("EVT", "Normal"),
                   nex = c(varex, varexnorm),
                   propex = c(varex/length(out_loss)*100, varexnorm/length(out_loss)*100))
colnames(dfex) <- c("Modelo", "Violações", "Proporção")
# stargazer(dfex, out = "artigo-apresentacao-tabela.tex", 
#           summary = FALSE, rownames = FALSE, font.size = "tiny",
#           style = "aer")
knitr::kable(dfex, format = "pandoc")


# Backtesting com refit ----------------------------------------------
# Primeiro um teste apenas para o Ibovespa
ibov.xts <- -assets.tbl$ts[[1]]
window.size <- ndays(ibov.xts[paste0("/", end)])
n.roll <- ndays(ibov.xts[paste0(backstart, "/")])
ibov_os <- (window.size+1):(window.size+n.roll)
realized <- ibov.xts[ibov_os] # Perdas realizadas durante o periodo fora da amostra

###### Modelo EVT condicional #########################################
rollspec <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                       variance.model = list(model = "eGARCH", garchOrder = c(2,1)),
                       distribution.model = "sstd")
# Ajusta os dados fora da amostra para o modelo EVT condicional
# Retorna um xts com os parametros da GPD e as medidas de risco
# CUIDADO!! uma rodada desta com 1000 observacoes fora da amostra
# pode levar mais de 6 HORAS (estimado 22 segundos para cada rolagem)
ibov_os.xts <- roll_fit_cevt(ibov.xts, rollspec, n.roll, window.size)
object.size(ibov_os.xts) # Retorna o tamanho em bytes, nao eh grande
# Salva os dados obtidos para tratamento posterior
saveRDS(ibov_os.xts, file = "ibov_os.xts.rds")
ibov_os.xts <- readRDS("ibov_os.xts.rds")
# Limpa a memoria para nao acumular muitos dados
#rm(list = "ibov_os.xts")

###### Modelo Normal incondicional #########################################
# Ajusta os dados para um modelo Normal incondicional
ibov_os_norm.xts <- roll_fit_unorm(ibov.xts, n.roll, window.size)
###### Modelo t-Student incondicional ######################################
# Ajusta os dados para um modelo t-Student incondicional
ibov_os_t.xts <- roll_fit_ut(ibov.xts, n.roll, window.size)
###### Modelo EVT incondicional #########################################
# Ajusta os dados para um modelo EVT incondicional





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

