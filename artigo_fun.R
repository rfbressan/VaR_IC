## Funcoes utilizadas em artigo.R
## artigo_fun.R
if(!all(c("rugarch", "fExtremes", "xts") %in% loadedNamespaces())){
  library(rugarch)
  library(fExtremes)
  library(xts)
} # Carrega os pacotes necessarios se faltantes

roll_fit <- function(data, spec, n.roll, window.size) {
  # Deve-se fazer o fit
  # Extrair residuos padronizados
  # Aplicar o gpfFit
  # Calcular zq e sq
  # Para cada interacao. n.roll vezes
  # Devolver tudo em um data.frame (ou lista)
  # data deve ser o xts com as perdas de todo o periodo
  tic <- Sys.time()
  NT <- length(data) # Tamanho total da amostra de dados
  cluster <- makePSOCKcluster(4)
  clusterEvalQ(cl = cluster, library(rugarch))
  clusterEvalQ(cl = cluster, library(fExtremes))
  clusterEvalQ(cl = cluster, library(xts))
  clusterExport(cluster, c("data", "spec", "n.roll", "window.size", "NT"),
                envir = environment())
  
  tmp.list <- parLapply(cl = cluster, 1:n.roll, fun = function(i){
    garch.fit <- ugarchfit(spec, 
                           data[(NT-(n.roll-i)-(window.size-1)):(NT-(n.roll-i))], 
                           solver = "hybrid")
    # Retirar os residuos padronizados, e o ultimo mu_t e sigma_t
    resid_z <- coredata(residuals(garch.fit, standardize = TRUE))
    mu_t <- last(fitted(garch.fit))
    sigma_t <- last(sigma(garch.fit))
    # Ajustar uma gpd aos dados de residuos padronizados
    gpd.fit <- gpdFit(resid_z, u = quantile(resid_z, 0.95))
    xi = gpd.fit@fit$par.ests[1]
    beta = gpd.fit@fit$par.ests[2]
    xi_se = gpd.fit@fit$par.ses[1]
    beta_se = gpd.fit@fit$par.ses[2]
    # Obter as medias de risco do residuo padronizado
    zq975 <- gpdRiskMeasures(gpd.fit, prob = 0.975)$quantile
    zq990 <- gpdRiskMeasures(gpd.fit, prob = 0.990)$quantile
    sq975 <- gpdRiskMeasures(gpd.fit, prob = 0.975)$shortfall
    sq990 <- gpdRiskMeasures(gpd.fit, prob = 0.990)$shortfall
    # Agora escalar e deslocar as medidas de risco de acordo com o modelo garch
    # Zq = mu_t+sigma_t*zq
    # Sq = mu_t+sigma_t*sq
    Zq975 <- mu_t+sigma_t*zq975
    Zq990 <- mu_t+sigma_t*zq990
    Sq975 <- mu_t+sigma_t*sq975
    Sq990 <- mu_t+sigma_t*sq990
    
    return(cbind(xi, beta, xi_se, beta_se, Zq975, Zq990, Sq975, Sq990))
  }) # fim do parLapply
  stopCluster(cluster)
  toc <- Sys.time()
  print(toc-tic)
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao 
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo em um xts indexado pelos dias fora da amostra
  ans <- xts(do.call(rbind, tmp.list),
             order.by = index(data[(NT-(n.roll-1)): NT]))
  colnames(ans) <- c("xi", "beta", "xi_se", "beta_se", "Zq975", "Zq990", "Sq975", "Sq990")
  return(ans)
  # Os valores retornados das medidas de risco devem ser comparadas
  # com os valores realizados NO DIA SEGUINTE a data onde foram calculadas
} # fim da roll_fit


