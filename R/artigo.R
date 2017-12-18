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
packages <- c("fExtremes", "rugarch", "xts", "PerformanceAnalytics", "xtable", "tidyverse", 
                      "broom", "purrr", "gridExtra", "ggplot2", "WeightedPortTest")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#library(fGarch)
#library(fBasics)
#library(QRM)
#library(timeSeries)
library(tidyverse)
library(broom)
library(purrr)
library(gridExtra)
library(ggplot2)
library(xtable)
library(WeightedPortTest)
#library(CADFtest) # Teste Dickey-Fuller
library(xts)
library(PerformanceAnalytics)
library(fExtremes)
library(rugarch)
source("./R/artigo_fun.R") # Carrega a funcao roll_fit para fazer o backtest

# AMOSTRA COM DADOS A PARTIR DE 31-08-2000
# 
start <- as.Date("2005-08-31")
end <- as.Date("2014-08-31")
backstart <- as.Date("2014-09-01")

list.returns <- function(asset, start) {
  tb <- read_csv(paste0("./input/artigo-", asset, ".csv"), 
                 col_types = cols_only(Date = col_date(), `Adj Close` = col_double()))
  prices <- xts(tb$`Adj Close`, order.by = tb$Date)[paste0(start, "/")]
  colnames(prices) <- "close"
  return(na.omit(Return.calculate(prices, method = "log")))
}
# Gera um tible com uma coluna com o codigo do ativo, a serie de retornos - ts e
# o nome do indice - id_name
assets <- c("BVSP", "GSPC", "GSPTSE", "IPSA", "MERV", "MXX")
lista <- lapply(assets, list.returns, start)
names(lista) <- assets
assets.tbl <- enframe(lista) %>% 
  bind_cols(tibble(id_name = c("IBovespa", "S&P500", "S&P TSE", "IPSA", "Merval", "IPC")))

colnames(assets.tbl) <- c("indice", "ts", "id_name")
# Qual o tamanho da janela in sample?
insample <-  enframe(map_dbl(lista, ~ndays(.x[paste0("/", end)])))
names(insample) <- c("indice", "insample")
# Remove a variavel lista e assets que agora sao desnecessarios
rm(lista)

# Estatisticas descritivas retornos-----------------------------------------
df.descritivas <- assets.tbl %>% 
  transmute(media = map_dbl(ts, ~ mean(.x)),
            mediana = map_dbl(ts, ~median(.x)),
            maximo = map_dbl(ts, ~max(.x)),
            minimo = map_dbl(ts, ~min(.x)),
            desvp = map_dbl(ts, ~sd(.x)),
            assim = map_dbl(ts, ~skewness(.x)),
            curtose = map_dbl(ts, ~kurtosis(.x)),
            jbstat = map_dbl(ts, ~jarqueberaTest(as.timeSeries(.x))@test$statistic),
            jbpvalue = map_dbl(ts, ~jarqueberaTest(as.timeSeries(.x))@test$p.value),
            q10stat = map_dbl(ts, ~Weighted.Box.test(.x, lag = 10, type = "Ljung-Box")$statistic),
            q10pvalue = map_dbl(ts, ~Weighted.Box.test(.x, lag = 10, type = "Ljung-Box")$p.value),
            q2_10stat = map_dbl(ts, ~Weighted.Box.test(.x, lag = 10, type = "Ljung-Box", sqrd.res = TRUE)$statistic),
            q2_10pvalue = map_dbl(ts, ~Weighted.Box.test(.x, lag = 10, type = "Ljung-Box", sqrd.res = TRUE)$p.value),
            nobs = map_int(ts, ~as.integer(length(.x)))) %>% 
  t() %>% 
  as.data.frame()
param <- c("Média", "Mediana", "Máximo", "Mínimo", "Desvp", "Assimetria", "Curtose",
                               "Jarque-Bera", "", "Q(10)", "", "$Q^2(10)$", "", "N.obs")
df.descritivas <- cbind(param, df.descritivas)
colnames(df.descritivas) <- c("Descritivas", assets.tbl$id_name)

# Cria o xtable
tab1 <- xtable(df.descritivas, caption = "Estatísticas descritivas dos retornos 
               (amostra completa de 31/08/2003 a 30/08/2017).",
               digits = 5,
               label = "tab:descritivas",
               auto = TRUE)
print.xtable(tab1, 
             file = "./tables/artigo-tab-descritivas.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

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
jpeg(filename = "./figs/artigo-qqplots.jpeg",
     width = 600, height = 800, quality = 100)
op <- par(mfrow = c(3,2),
          mar = c(4, 3, 3, 2))
for(i in 1:dim(assets.tbl)[1]){
  qqnormPlot(assets.tbl$ts[[i]], labels = FALSE, title = FALSE, mtext = FALSE,
             main = assets.tbl$id_name[[i]], 
             xlab = "")
             #ylab = "Amostra")
}
par(op)
dev.off()

# Modelo eGARCH in Sample-----------------------------------------------------------
## Ja as distribuicoes de zt a normal e t-Student nao apresentam bom fit
## A Johnson, GED, NIG, SkewStudent e a Ghyp sao melhores
## Lembrando, o modelo das perdas é AR(1) e a volatilidade é eGARCH(2,1)
## L_t=mu_t+e_t      mu_t=mu+ar1*mu_t-1+e_t
## e_t=sigma_t*z_t   ln(sigma^2_t)=omega+alpha1*z_t-1+gamma1(|z_t-1|-E[|zt-1|])+beta1*ln(sigma^2_t-1)
## LEMBRAR: ts contem os retornos e não as perdas!

# Uma especificacao para cada ativo
ruspec <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "eGARCH", garchOrder = c(2,1)),
                     distribution.model = "norm")
garch.specs <- replicate(length(assets), ruspec)
names(garch.specs) <- assets
garch.specs <- enframe(garch.specs)
colnames(garch.specs) <- c("indice", "spec")
## Modelando as PERDAS!! parametro eh loss para a funcao ugarchfit
# Deixamos um numero outsample para fazer o backtesting. Diferente para cada ativo
# garch.models vai conter o modelo eGARCH dentro da amostra. Apresentar os 
# parametros e seus erros padrao robustos
# Depois apresentar novamente estatisticas como JB, Q e Q^2 para os erros padronizados
garch.models <- assets.tbl[,1:3] %>% 
  inner_join(insample, by = "indice") %>% 
  inner_join(garch.specs, by = "indice") %>% 
  mutate(loss = map(ts, ~-.x),
         loss_in = map2(ts, insample, ~-.x[1:.y]),
         garch_fit = map2(spec, loss_in, ~ugarchfit(.x, 
                                                    .y,
                                                    solver = "hybrid")),
         ts = NULL)

# Mostra a convergencia para cada modelo
lapply(garch.models$garch_fit, convergence)

# Sumarios dos modelos GARCH
show(garch.models$garch_fit[[1]])

## Construindo a tabela com os parametros estimados do eGARCH in sample
garch.models.par <- garch.models %>% 
  transmute(mu = map(.$garch_fit, ~.x@fit$robust.matcoef["mu", c(1, 4)]), # Estimativa e P-valor
            ar1 = map(.$garch_fit, ~.x@fit$robust.matcoef["ar1", c(1, 4)]),
            omega = map(.$garch_fit, ~.x@fit$robust.matcoef["omega", c(1, 4)]),
            alpha1 = map(.$garch_fit, ~.x@fit$robust.matcoef["alpha1", c(1, 4)]),
            alpha2 = map(.$garch_fit, ~.x@fit$robust.matcoef["alpha2", c(1, 4)]),
            beta1 = map(.$garch_fit, ~.x@fit$robust.matcoef["beta1", c(1, 4)]),
            gamma1 = map(.$garch_fit, ~.x@fit$robust.matcoef["gamma1", c(1, 4)]),
            gamma2 = map(.$garch_fit, ~.x@fit$robust.matcoef["gamma2", c(1, 4)]),
            skew = map(.$garch_fit, ~.x@fit$robust.matcoef["skew", c(1, 4)]),
            shape = map(.$garch_fit, ~.x@fit$robust.matcoef["shape", c(1, 4)])) %>% 
  t() %>% 
  as.data.frame() %>% 
  unnest() 
  colnames(garch.models.par) <- garch.models$id_name
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
  garchcoef <- cbind(par = param, garch.models.par)
  colnames(garchcoef)[1] <- "Parâmetros"
  rm(garch.models.par) # Nao utilizado mais
  # Xtable
  tab2 <- xtable(garchcoef, caption = "Par\\^ametros estimados do modelo eGARCH. Valores p apresentados 
               de acordo com erros padrão robustos. (amostra de trabalho entre 31/08/2003 a 31/08/2013).",
                 digits = 5,
                 label = "tab:garchcoef",
                 auto = TRUE)
  print.xtable(tab2, 
               file = "./tables/artigo-tab-garchcoef.tex",
               caption.placement = "top",
               table.placement = "H",
               sanitize.colnames.function = NULL,
               sanitize.text.function = function(x) {x},
               include.rownames = FALSE)
##################### NAO USADO ###################################
## Para obter os valores de mu_t, sigma_t e z_t fora da amostra
## se utiliza o metodo ugarchfilter com os parametros fixados nos 
## valores estimados pelo ugarchfit
# garch.filtered <- garch.models %>% 
#   mutate(n_old = map_int(garch_fit, ~.x@model$modeldata$T),
#          f_params = map(garch_fit, ~as.list(coef(.x))),
#          out_spec = map(f_params, ~ugarchspec(mean.model = list(armaOrder = c(1,0)),
#                                               variance.model = list(model = "eGARCH", garchOrder = c(2,1)),
#                                               distribution.model = "sstd",
#                                               fixed.pars = .x))) %>% 
#   mutate(garch_filter = pmap(., ~ugarchfilter(..8,
#                                               ..4,
#                                               n.old = ..6)))
# 
# # Teste para verificar os residuos padronizados sao os mesmos inicialmente
# head(rugarch::residuals(garch.filtered$garch_fit[[2]], standardize = TRUE))
# head(rugarch::residuals(garch.filtered$garch_filter[[2]], standardize = TRUE))
# Agora garch.filtered contem tudo sobre os modelos Garch, nao precisamos mais de garch.models
##########################################################################################

# Terminamos ainda com efeito alavancagem a ser modelado para alguns casos
# signbias mostra significancia do efeito negativo
# Talvez um modelo GJR ou APARCH possa resolver. Nao resolveram, melhor eGARCH(2,1)

# Gerar 6 figuras com estes 4 graficos ACF
for(i in 1:dim(garch.models)[1]) {
  jpeg(filename = paste0("./figs/artigo-acf-", garch.models$id_name[i], ".jpeg"),
       width = 800, height = 800, quality = 100)
  op <- par(mfrow=c(2,2))
  plot(garch.models$garch_fit[[i]], which = 4)
  plot(garch.models$garch_fit[[i]], which = 5)
  plot(garch.models$garch_fit[[i]], which = 10)
  plot(garch.models$garch_fit[[i]], which = 11)
  par(op)
  dev.off()
}
file.rename(c("./figs/artigo-acf-S&P500.jpeg", "./figs/artigo-acf-S&P TSE.jpeg"), 
            c("./figs/artigo-acf-SP500.jpeg", "./figs/artigo-acf-SP-TSE.jpeg"))
## Estatisticas modelo eGARCH in sample
# JB, Q e Q^2 para os residuos padronizados
garch.models.stats <- garch.models %>%
  transmute(resid_z = map(.$garch_fit, ~residuals(.x, standardize = TRUE)),
            jbstat = map_dbl(resid_z, ~jarqueberaTest(as.timeSeries(.x))@test$statistic),
            jbpvalue = map_dbl(resid_z, ~jarqueberaTest(as.timeSeries(.x))@test$p.value),
            q10stat = map_dbl(resid_z, ~Weighted.Box.test(.x, lag = 10, 
                                                      type = "Ljung-Box")$statistic),
            q10pvalue = map_dbl(resid_z, ~Weighted.Box.test(.x, lag = 10, 
                                                          type = "Ljung-Box")$p.value),
            q2_10stat = map_dbl(resid_z, ~Weighted.Box.test(.x, lag = 10, 
                                                            type = "Ljung-Box", 
                                                            sqrd.res = TRUE)$statistic),
            q2_10pvalue = map_dbl(resid_z, ~Weighted.Box.test(.x, lag = 10, 
                                                            type = "Ljung-Box", 
                                                            sqrd.res = TRUE)$p.value)) %>% 
  mutate(resid_z = NULL) %>% 
  t() %>% 
  as.data.frame()
param <- c("Jarque-Bera", "",
           "Q(10)", "",
           "$Q^2(10)$", "")
garch.models.stats <- cbind(param, garch.models.stats)
colnames(garch.models.stats) <- c("Estatística", garch.models$id_name)
# Xtable
tab3 <- xtable(garch.models.stats, caption = "Estatísticas de diagnóstico para o modelo eGARCH. 
               (amostra de trabalho entre 31/08/2005 a 31/08/2013).",
               digits = 5,
               label = "tab:garchstats",
               auto = TRUE)
print.xtable(tab3, 
             file = "./tables/artigo-tab-garchstats.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

# Modelo EVT para os residuos padronizados in Sample --------------------------------

## Os residuos podem ser retirados atraves do metodo residuals com a opcao "standardize=T"
## os valores zq (quantile) e sq (shortfall)
evt.models <- garch.models %>% 
  transmute(indice = indice,
            id_name = id_name,
            loss_in = loss_in,
            insample = insample,
            resid_z = map(garch_fit, ~coredata(residuals(.x, standardize = TRUE))),
            mut = map(garch_fit, ~coredata(fitted(.x))),
            sigmat = map(garch_fit, ~coredata(sigma(.x))),
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
testez <- coredata(residuals(garch.models$garch_fit[[1]], standardize = TRUE))
teste_evir <- gpd(testez, threshold = quantile(testez, 0.925)) # Mesmos valores do fExtremes
meplot(testez, type = "l")
# shape(testez, models = 10, start = 90, end = 150)
# ## Teste com o pacote evd
# teste_evd <- fpot(testez, quantile(testez, 0.95))
# fitted.values(teste_evd)
# std.errors(teste_evd)  # Iguais ao fExtremes e evir
# mrlplot(testez)
evtcoef <- evt.models %>% 
  transmute(insample = insample,
            u = u,
            Nu = Nu,
            xi = xi,
            xi_se = xi_se,
            beta = beta,
            beta_se = beta_se,
            zq975 = zq975,
            zq990 = zq990) %>% 
  t() %>% 
  as.data.frame()
param <- c("Obs. dentro amostra", "Limiar", "Número de excessos", "Parâmetro forma GPD", "Erro padrão",
           "Parâmetro escala GPD", "Erro padrão", "Quantil 97.5\\%", "Quantil 99.0\\%")
evtcoef <- cbind(param, evtcoef)
colnames(evtcoef) <- c("", evt.models$id_name)
# Xtable
tab4 <- xtable(evtcoef, caption = "Parâmetros estimados para o modelo EVT dos resíduos padronizados. 
               (amostra de trabalho entre 31/08/2003 a 31/08/2013).",
               digits = 5,
               label = "tab:evtcoef",
               auto = TRUE)
print.xtable(tab4, 
             file = "./tables/artigo-tab-evtcoef.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

## Gráficos para analisar a qualidade do gpdFit
op <- par(mfrow=c(2,2))
plot(evt.models$gpdfit[[1]], which='all')
par(op)

# Reconstruindo o VaR e o ES condicionais in Sample ---------------------------------

## VaR: xq_t = mu_t+1 + sigma_t+1*zq
## ES: Sq_t = mu_t+1 + sigma_t+1*sq
riskmeasures <- evt.models %>%
  transmute(indice = indice,
            id_name = id_name,
            loss_in = loss_in,
            VaR975 = pmap(., ~(..6+..7*..15)), ## VaR = mu_t+1 + sigma_t+1*zq
            VaR990 = pmap(., ~(..6+..7*..16)),
            ES975 = pmap(., ~(..6+..7*..17)),  ## ES = mu_t+1 + sigma_t+1*sq
            ES990 = pmap(., ~(..6+..7*..18)))
# %>% 
#   mutate(out_VaR975 = pmap(., ~..6[c((..5+1):(..5+..4-1))]), #out_VaR = VaR[(n_old+1):(n_old+out-1)]
#          out_VaR990 = pmap(., ~..7[c((..5+1):(..5+..4-1))]),
#          out_ES975 = pmap(., ~..8[c((..5+1):(..5+..4-1))]),  #out_ES = ES[(n_old+1):(n_old+out-1)]
#          out_ES990 = pmap(., ~..9[c((..5+1):(..5+..4-1))]),
#          out_loss = map2(loss, n_old, ~coredata(.x[-c(1:(.y+1))])[, 1, drop = TRUE])) #out_los = loss[-c(1:n_old+1)]
# Plotando os valores dentro da amostra
plot_risks <- function(loss, VaR, ES, id_name) {
  xts <- merge(loss = loss, VaR = VaR, ES = ES)
  plot <- ggplot(xts, aes(x = Index))+
    geom_line(aes(y = loss), color = "black")+
    geom_line(aes(y = VaR), color = "red")+
    geom_line(aes(y = ES), color = "darkgreen")+
    labs(x = "", y = "", title = id_name)
  return(plot)
}
VaR_plots <- riskmeasures %>% 
  transmute(VaR975_plot = pmap(., ~plot_risks(..3, ..4, ..6, ..2)),
            VaR990_plot = pmap(., ~plot_risks(..3, ..5, ..7, ..2)))
VaR_plots$VaR975_plot[[1]]+
  coord_cartesian(xlim = c(as.Date("2011-08-31"), 
                           as.Date("2014-08-31")))
# Verifica quantas violacoes
sum(riskmeasures$loss_in[[1]] > riskmeasures$VaR975[[1]])

grid.arrange(grobs = VaR_plots$VaR975_plot)
grid.arrange(grobs = VaR_plots$VaR990_plot)

# Backtesting com refit ----------------------------------------------
# Monta o tibble para as estimacoes out of sample
assets_os.tbl <- garch.models %>% 
  transmute(indice = indice,
            id_name = id_name,
            loss = loss,
            window.size = map_int(loss, ~as.integer(ndays(.x[paste0("/", end)]))),
            n.roll = map2_int(loss, window.size, ~length(.x)-.y),
            spec = spec)
# Periodo fora da amostra + 1, para fazer as comparacoes entre
# VaRt e realizado t+1
# realized tera n.roll-1 observacoes devido ao deslocamento das observacoes de VaR e ES
# para serem comparados.
# As series de VaR e ES no tibble de risco tambem terao n.roll-1 observacoes e estarao nas 
# datas corretas para COMPARACAO (a medida foi calculada no dia anterior)
realized <- assets_os.tbl %>% 
  transmute(indice = indice,
            real = pmap(., ~..3[(..4+2):(..4+..5)])) # real = loss[(window.size+2):(window.size+n.roll)]

######################### NAO USADO ####################################################################         
## Modelo EVT condicional ###
#rollspec <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
#                       variance.model = list(model = "eGARCH", garchOrder = c(2,1)),
#                       distribution.model = "norm")
# Ajusta os dados fora da amostra para o modelo EVT condicional
# Retorna um xts com os parametros da GPD e as medidas de risco
# CUIDADO!! uma rodada desta com 1000 observacoes fora da amostra
# pode levar mais de 6 HORAS (estimado 22 segundos para cada rolagem)
# ibov_os_cevt.list <- roll_fit_cevt(ibov.xts, window.size, n.roll, rollspec)
# object.size(ibov_os.list) # Retorna o tamanho em bytes, nao eh grande
# Salva os dados obtidos para tratamento posterior
# saveRDS(ibov_os.list, file = "ibov_os.list.rds")
# ibov_os.xts <- readRDS("ibov_os.list.rds")
# Limpa a memoria para nao acumular muitos dados
#rm(list = "ibov_os.list")

### Modelo EVT incondicional ###
# Ajusta os dados para um modelo EVT incondicional
#ibov_os_uevt.xts <- roll_fit_uevt(ibov.xts, window.size, n.roll, rollspec)
###### Modelo Normal incondicional ##
# Ajusta os dados para um modelo Normal incondicional
#ibov_os_norm.xts <- roll_fit_unorm(ibov.xts, window.size, n.roll)
###### Modelo t-Student incondicional ###
# Ajusta os dados para um modelo t-Student incondicional
#ibov_os_t.xts <- roll_fit_ut(ibov.xts, window.size, n.roll)
###### Modelo RiskMetrics com suavizacao exponencial ###
# Ajusta os dados para um modelo RiskMetrics
#ibov_os_riskmetrics <- roll_fit_riskmetrics(ibov.xts, window.size, n.roll)
##########################################################################################

# Testes estatisticos para o VaR ------------------------------------------
# Quais testes fazer?
# VaRTest possui 2 testes, incondicional de Kupiec1995 e condicional de Christoffersen2001
# que nao esta convergindo
# Existe tambem o teste de Christoffersen2004 de tempo entre as violacoes
# Pode ser feito um teste do tipo Model Confidence Set - MCS - com a funcao VaRloss e 
# mcsTest

######### TESTE PARA OS 6 INDICES ###
# roll_fit(data, window.size, n.roll, spec, models)
# de assets_os.tbl temos os seguintes numeros
# 1: indice 
# 2: id_name
# 3: loss 
# 4: window.size 
# 5: n.roll 
# 6: spec 
models <- c("cevt", "unorm", "ut", "uevt", "riskmetrics")
#teste.tbl <- assets_os.tbl[2,] # Copia os dados para um tibble de teste
## ATENCAO! Aqui eh alterado o valor de n.roll para o teste ser rapido
## no artigo completo deve-se EXCLUIR ESTA PROXIMA LINHA
#teste.tbl$n.roll <- 100
#teste.tbl$n.roll <- as.integer(rep(15, 6))
cat("\nInicio do map roll_fit:", as.character(Sys.time()))
os_roll.tbl <- assets_os.tbl %>% 
  transmute(indice = indice,
            id_name = id_name,
            roll.fit = pmap(., ~roll_fit(..3, ..4, ..5, ..6, models)))
cat("\nFim do map roll_fit:", as.character(Sys.time()))
saveRDS(os_roll.tbl, file = "os_roll_tbl.rds")

os_roll_unnest <- os_roll.tbl %>% 
  unnest() %>% 
  mutate(risk.tbl = map(roll, ~.x$risk.tbl),
         param.xts = map(roll, ~.x$param.xts))

# Extrai a evolucao dos parametros de cada modelo no tempo
param.tbl <- os_roll_unnest %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            param.xts = param.xts)
# Extrai a evolucao das medidas de risco de cada modelo, para cada cobertura, no tempo
os_risk.tbl <- os_roll_unnest %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            risk.tbl = risk.tbl) %>% 
  unnest()
format(object.size(os_risk.tbl), units = "Kb") # Verifica o tamanho do objeto

## Faz os testes de VaR para os 6 indices
## VaRTest(alpha = 0.05, actual, VaR, conf.level = 0.95)
# teste_realized <- realized %>% # Copia os dados para um tibble de teste
#   inner_join(teste.tbl, by = "indice") %>% 
#   transmute(indice = indice,
#             real = map2(real, n.roll, ~.x[1:(.y-1)]))
vartest.tbl <- os_risk.tbl %>% 
  left_join(realized, by = "indice") %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            VaRtest = pmap(., ~VaRTest(..4, -coredata(..7), -coredata(..5)))) # Troca o sinal!!

# Plot do VaR e violacoes
lin <- 2 # Escolhe uma das linhas de risk, de 1 a 10 para ficar no Ibovespa
VaRplot(os_risk.tbl$coverage[lin], -realized$real[[1]], -os_risk.tbl$VaR.xts[[lin]])

## VaRDurTest
vardurtest.tbl <- os_risk.tbl %>% 
  left_join(os_realized, by = "indice") %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            VaRDurTest = pmap(., ~VaRDurTest(..4, -coredata(..7), -coredata(..5)))) # Troca o sinal!!

## MCS para os VaR atraves da funcao VaRLoss
mcs.tbl <- os_risk.tbl %>% 
  left_join(realized, by = "indice") %>% 
  transmute(indice = indice,
            model_type = model_type,
            coverage = coverage,
            VaRloss = pmap(., ~VaRloss(..4, -coredata(..7), -coredata(..5)))) %>% # Troca o sinal!!
  group_by(indice, coverage) %>% 
  summarise(loss_matrix = list(do.call(cbind, VaRloss))) %>% 
  mutate(mcs_test = map2(loss_matrix, coverage, ~mcsTest(.x, .y)))

# Testes estatisticos para o ES -------------------------------------------
# ATENCAO!! Este teste so faz sentido ser realizados apos os teste de VaR
# aprovarem o modelo
# ESTest(alpha = 0.05, actual, ES, VaR, conf.level = 0.95, boot = FALSE, n.boot = 1000)
# Pelo codigo fonte nao parece ser bem o teste implementado pelo McNeil2000 com o algoritmo
# do Efron1993
ru_estest.tbl <- os_risk.tbl %>% 
  left_join(realized, by = "indice") %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            EStest = pmap(., ~ESTest(..4, -coredata(..7), -coredata(..6), -coredata(..5), # troca o sinal!!
                                     boot = TRUE, n.boot = 1000)))

# ATENCAO!! Este teste so faz sentido ser realizados apos os teste de VaR
# aprovarem o modelo
# Chamando a funcao propria es_test
es.test <- teste_risk.tbl %>%
  left_join(teste_realized, by = "indice") %>%
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            EStest = pmap(., ~es_test(..4, coredata(..7), coredata(..6), coredata(..5)))) # Nao troca o sinal

losses <- coredata(riskmeasures$loss_in[[2]])
VaR <- riskmeasures$VaR975[[2]]
ES <- riskmeasures$ES975[[2]]
t_es.test <- es_test(0.0275, losses, ES, VaR, n.boot = 100)


