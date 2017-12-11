## Funcoes utilizadas em artigo.R
## artigo_fun.R
if(!all(c("rugarch", "fExtremes", "xts") %in% loadedNamespaces())){
  library(rugarch)
  library(fExtremes)
  library(xts)
} # Carrega os pacotes necessarios se faltantes

cores <- detectCores() # Quantos cores estao rodando

# roll_fit_cevt -----------------------------------------------------------
roll_fit_cevt <- function(data, spec, n.roll, window.size) {
  # Deve-se fazer o fit
  # Extrair residuos padronizados
  # Aplicar o gpfFit
  # Calcular zq e sq
  # Para cada interacao. n.roll vezes
  # Devolver tudo em um data.frame
  # data deve ser o xts com as perdas de todo o periodo
  tic <- Sys.time()
  NT <- length(data) # Tamanho total da amostra de dados
  # cluster <- makePSOCKcluster(4)
  # clusterEvalQ(cl = cluster, library(rugarch))
  # clusterEvalQ(cl = cluster, library(fExtremes))
  # clusterEvalQ(cl = cluster, library(xts))
  # clusterExport(cluster, c("data", "spec", "n.roll", "window.size", "NT"),
  #               envir = environment())

  tmp.list <- mclapply(1:n.roll, function(i){
    garch.fit <- try(ugarchfit(spec,
                           data[i:(window.size+i)],
                           solver = "hybrid"))
    # 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
    if(inherits(garch.fit, 'try-error') || convergence(garch.fit)!=0 || is.null(garch.fit@fit$cvar)){
      lapply_ans <- t(cbind(rep(NA, 9)))
      # Se algum erro, envia como resposta NA, mas nao interrompe a estimacao
      # Depois no xts retornado pela funcao roll_fit, preencher os NA com os dados
      # da estimacao anterior com na.locf
    } else{
      # Salvar o coeficiente beta1 apenas para mostrar a evolucao deste ao longo
      # do período
      beta1 <- coef(garch.fit)["beta1"]
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

      # Verbose mode
      #cat("Estimacao numero:", i, "as", as.character(Sys.time()))
      lapply_ans <- cbind(beta1, xi, beta, xi_se, beta_se, Zq975, Zq990, Sq975, Sq990)
    } # fim do else
    return(lapply_ans)
  },
  mc.cores = cores) # fim do mclapply
  #stopCluster(cluster)
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]
  ans <- xts(ans, order.by = index(data[(window.size+2):(window.size+n.roll)]))
  colnames(ans) <- c("beta1", "xi", "beta", "xi_se", "beta_se", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  toc <- Sys.time()
  print(toc-tic)
  return(ans)
  # Os valores retornados das medidas de risco devem ser comparadas
  # com os valores realizados NO DIA SEGUINTE a data onde foram calculadas
} # fim da roll_fit_cevt


# roll_fit_unorm ----------------------------------------------------------
# Ajusta os dados para um modelo Normal incondicional
roll_fit_unorm <- function(data, n.roll, window.size) {
  tic <- Sys.time()
  tmp.list <- mclapply(1:n.roll, function(i){
    xts <- data[i:(window.size+i)]
    mean <- mean(xts)
    sd <- sd(xts)
    Zq975 <- qnorm(0.975, mean, sd)
    Zq990 <- qnorm(0.990, mean, sd)
    # Equacao do ES retirada de Pfaff2013 p. 36, eq. 4.5
    # ES_a = 1/(1-a) * int_a^1 q_l(x)dx
    Sq975 <- integrate(function(x){
      qnorm(x, mean = mean, sd = sd)},
      0.975,
      1)$value / (1-0.975)
    Sq990 <- integrate(function(x){
      qnorm(x, mean = mean, sd = sd)},
      0.990,
      1)$value / (1-0.990)
    return(cbind(mean, sd, Zq975, Zq990, Sq975, Sq990))
  },
  mc.cores = cores) # Fim do lapply
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao 
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem 
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]           
  ans <- xts(ans, order.by = index(data[(window.size+2):(window.size+n.roll)]))
  colnames(ans) <- c("mu", "sigma", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  toc <- Sys.time()
  print(toc-tic)
  return(ans)
}

# roll_fit_ut ----------------------------------------------------------
# Ajusta os dados para um modelo t-Student incondicional
roll_fit_ut <- function(data, n.roll, window.size) {
  tic <- Sys.time()
  tmp.list <- mclapply(1:n.roll, function(i){
    xts <- data[i:(window.size+i)]
    t_fit <- fitdist(distribution = "std", xts)
    t_mu <- t_fit$pars["mu"]
    t_sigma <- t_fit$pars["sigma"]
    t_shape <- t_fit$pars["shape"]
    Zq975 <- qdist(distribution = "std", p = 0.975,
                   mu = t_mu,
                   sigma = t_sigma,
                   shape = t_shape)
    Zq990 <- qdist(distribution = "std", p = 0.990,
                   mu = t_mu,
                   sigma = t_sigma,
                   shape = t_shape)
    # Equacao do ES retirada de Pfaff2013 p. 36, eq. 4.5
    # ES_a = 1/(1-a) * int_a^1 q_l(x)dx
    Sq975 <- integrate(function(x){
      qdist(distribution = "std", x, mu = t_mu, sigma = t_sigma, shape = t_shape)},
      0.975,
      1)$value / (1-0.975)
    Sq990 <- integrate(function(x){
      qdist(distribution = "std", x, mu = t_mu, sigma = t_sigma, shape = t_shape)},
      0.990,
      1)$value / (1-0.990)
    return(cbind(t_mu, t_sigma, t_shape, Zq975, Zq990, Sq975, Sq990))
  },
  mc.cores = cores) # Fim do mclapply
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao 
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem 
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]           
  ans <- xts(ans, order.by = index(data[(window.size+2):(window.size+n.roll)]))
  colnames(ans) <- c("mu", "sigma", "shape", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  toc <- Sys.time()
  print(toc-tic)
  return(ans)
}

# roll_fit_uevt -----------------------------------------------------------
roll_fit_uevt <- function(data, spec, n.roll, window.size) {
  # Deve-se fazer o fit
  # Extrair residuos padronizados
  # Aplicar o gpfFit
  # Calcular zq e sq com base na media e variancia INCONDICIONAIS
  # Para cada interacao. n.roll vezes
  # Devolver tudo em um data.frame (ou lista)
  # data deve ser o xts com as perdas de todo o periodo
  tic <- Sys.time()
  NT <- length(data) # Tamanho total da amostra de dados
  # cluster <- makePSOCKcluster(4)
  # clusterEvalQ(cl = cluster, library(rugarch))
  # clusterEvalQ(cl = cluster, library(fExtremes))
  # clusterEvalQ(cl = cluster, library(xts))
  # clusterExport(cluster, c("data", "spec", "n.roll", "window.size", "NT"),
  #               envir = environment())
  # 
  tmp.list <- mclapply(1:n.roll, function(i){
    garch.fit <- try(ugarchfit(spec, 
                               data[i:(window.size+i)], 
                               solver = "hybrid"))
    # 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
    if(inherits(garch.fit, 'try-error') || convergence(garch.fit)!=0 || is.null(garch.fit@fit$cvar)){
      lapply_ans <- t(cbind(rep(NA, 10)))
      # Se algum erro, envia como resposta NA, mas nao interrompe a estimacao
      # Depois no xts retornado pela funcao roll_fit, preencher os NA com os dados 
      # da estimacao anterior com na.locf
    } else{
      # Salvar o coeficiente beta1 apenas para mostrar a evolucao deste ao longo
      # do período
      beta1 <- coef(garch.fit)["beta1"]
      # Retirar os residuos padronizados, e o ultimo mu_t e sigma_t
      resid_z <- coredata(residuals(garch.fit, standardize = TRUE))
      mu_t <- uncmean(garch.fit)
      sigma_t <- uncvariance(garch.fit)
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
      
      # Verbose mode
      #cat("Estimacao numero:", i, "as", as.character(Sys.time()))
      lapply_ans <- cbind(beta1, xi, beta, xi_se, beta_se, Zq975, Zq990, Sq975, Sq990)
    } # fim do else
    return(lapply_ans)
  },
  mc.cores = cores) # fim do mclapply
  #stopCluster(cluster)
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao 
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem 
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]           
  ans <- xts(ans, order.by = index(data[(window.size+2):(window.size+n.roll)]))
  colnames(ans) <- c("beta1", "xi", "beta", "xi_se", "beta_se", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  toc <- Sys.time()
  print(toc-tic)
  return(ans)
  # Os valores retornados das medidas de risco devem ser comparadas
  # com os valores realizados NO DIA SEGUINTE a data onde foram calculadas
} # fim da roll_fit_cevt

