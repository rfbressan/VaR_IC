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
cap <- paste("Estatísticas descritivas dos retornos (amostra completa de",
           format(start+1, "%d/%m/%Y"), "a 30/08/2017).")
tab1 <- xtable(df.descritivas, 
               caption = cap,
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
# adf <- lapply(assets.tbl$ts, CADFtest)
# names(adf) <- assets.tbl$indice
# adf <- enframe(adf)
# adf <- adf %>%
#   mutate(teste = map(value, ~tidy(.x)))
# adf$value <- NULL
# adf <- unnest(adf)
# adf

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
            gamma2 = map(.$garch_fit, ~.x@fit$robust.matcoef["gamma2", c(1, 4)])) %>% 
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
             "$\\gamma_2$", "")
  garchcoef <- cbind(par = param, garch.models.par)
  colnames(garchcoef)[1] <- "Parâmetros"
  rm(garch.models.par) # Nao utilizado mais
  # Xtable
  cap <- paste("Par\\^ametros estimados do modelo eGARCH. Valores p apresentados de acordo 
com erros padrão robustos. (amostra de trabalho entre",
               format(start+1, "%d/%m/%Y"), "a",
               format(end, "%d/%m/%Y"), ").")
  tab2 <- xtable(garchcoef, 
                 caption = cap,
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
cap <- paste("Estatísticas de diagnóstico para o modelo eGARCH. 
               (amostra de trabalho entre",
             format(start+1, "%d/%m/%Y"), "a",
             format(end, "%d/%m/%Y"), ").")

tab3 <- xtable(garch.models.stats, 
               caption = cap,
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
cap <- paste("Parâmetros estimados para o modelo EVT dos resíduos padronizados. 
               (amostra de trabalho entre",
             format(start+1, "%d/%m/%Y"), "a",
             format(end, "%d/%m/%Y"), ").")

tab4 <- xtable(evtcoef, 
               caption = cap,
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
for(i in seq_len(dim(garch.models)[1])) {
  jpeg(filename = paste0("./figs/artigo-evtgof-", garch.models$indice[i], ".jpeg"),
       width = 800, height = 800, quality = 100)
  op <- par(mfrow=c(2,2))
  plot(evt.models$gpdfit[[i]], which='all')
  par(op)
  dev.off()
}

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
models <- c("cevt", "cnorm", "ct", "uevt", "unorm", "ut", "riskmetrics")
# teste_assets_os <- readRDS("./input/teste_assets_os.rds") # Copia dados de teste
# teste_realized <- readRDS("./input/teste_realized.rds")
# teste_os_roll <- teste_assets_os %>% 
#   transmute(indice = indice,
#             id_name = id_name,
#             roll.fit = pmap(., ~roll_fit(..3, ..4, ..5, ..6, models)))

## ATENCAO! Aqui eh alterado o valor de n.roll para o teste ser rapido
cat("\nInicio do map roll_fit:", as.character(Sys.time()))
# os_roll.tbl <- assets_os.tbl %>% 
#   transmute(indice = indice,
#             id_name = id_name,
#             roll.fit = pmap(., ~roll_fit(..3, ..4, ..5, ..6, models)))
# cat("\nFim do map roll_fit:", as.character(Sys.time()))
# saveRDS(os_roll.tbl, file = "./output/os_roll_tbl.rds")
os_roll.tbl <- readRDS(file = "./output/os_roll_tbl.rds")

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

## Tabela com os percentuais de violacoes
varviolations.tbl <- os_risk.tbl %>% 
  left_join(realized, by = "indice") %>% 
  mutate(violations = map2_dbl(VaR.xts, real, ~100*sum(.y > .x)/length(.y)),
         cov = 100*coverage) %>% 
  select(id_name, cov, model_type, violations) %>% 
  spread(key = id_name, value = violations) %>% 
  mutate(model_type = map_chr(model_type, ~switch(.x,
                                              cevt = "EVT Condicional",
                                              cnorm = "Normal Condicional",
                                              ct = "t-Student Condicional",
                                              uevt = "EVT Incond. Filtrada",
                                              unorm = "Normal Incondicional",
                                              ut = "t-Student Incondicional",
                                              riskmetrics = "RiskMetrics") # Fim do switch
  )) # Fim do map_chr e mutate
colnames(varviolations.tbl)[2] <- "Modelo"
varviolations.tbl <- add_row(varviolations.tbl, cov = 1.0, Modelo = "Cobertura = 1\\%", .before = 1)
varviolations.tbl <- add_row(varviolations.tbl, cov = 2.5, Modelo = "Cobertura = 2.5\\%", .before = 9)
varviolations.tbl$cov <- NULL # Retira a coluna cov, que nao eh mais necessaria

# Xtable
cap <- paste("Percentual de violações. (fora da amostra, dados entre",
             format(backstart+1, "%d/%m/%Y"), "e 31/08/2017")

tab5 <- xtable(varviolations.tbl, 
               caption = cap,
               digits = 2,
               label = "tab:varviol",
               auto = TRUE)
print.xtable(tab5, 
             file = "./tables/artigo-tab-varviol.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

## Faz os testes de VaR para os 6 indices
## VaRTest(alpha = 0.05, actual, VaR, conf.level = 0.95)
## var_test(cover = 0.025, loss, var, conf.level = 0.95)
vartest.tbl <- os_risk.tbl %>% 
  left_join(realized, by = "indice") %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            VaRtest = pmap(., ~var_test(..4, coredata(..7), coredata(..5)))) %>% 
  unnest() %>% 
  select(id_name, coverage, model_type, uc.LRstat, uc.LRp, cc.LRstat, cc.LRp) %>% 
  gather(key = stat_name, value = stat_value, -c(id_name, coverage, model_type), factor_key = TRUE) %>% 
  spread(key = id_name, value = stat_value)

###### Teste da funcao var_test
teste_risk <- os_risk.tbl %>% 
  subset(subset = (indice == "IPSA" & model_type == "cevt" & coverage == 0.01)) %>% 
  left_join(realized, by = "indice")
var_test.tbl <- var_test(cover = teste_risk$coverage,
                         loss = unlist(teste_risk$real), 
                         var = unlist(teste_risk$VaR.xts))

# Plot do VaR e violacoes
teste_varplot <- os_risk.tbl %>% 
  subset(subset = (indice == "IPSA" & model_type == "unorm" & coverage == 0.01)) %>% 
  left_join(realized, by = "indice")
VaRplot(teste_varplot$coverage, -teste_varplot$real[[1]], -teste_varplot$VaR.xts[[1]])

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


