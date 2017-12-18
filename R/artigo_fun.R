## Funcoes utilizadas em artigo.R
## artigo_fun.R
if(!all(c("rugarch", "fExtremes", "xts", "tidyverse") %in% loadedNamespaces())){
  library(tidyverse)
  library(xts)
  library(fExtremes)
  library(rugarch)
} # Carrega os pacotes necessarios se faltantes

cores <- detectCores() # Quantos cores estao rodando

# roll_fit ----------------------------------------------------------------

roll_fit <- function(data, window.size, n.roll, spec, models) {
  cat("\nroll_fit do ativo inicio:", as.character(Sys.time()))
  # Check para os argumentos
  if(!is.xts(data)) stop("roll_fit: data deve ser um xts")
  if(any(is.null(window.size), is.null(n.roll), is.null(spec), is.null(models)))
    stop("roll_fit: Devem ser passados todos os argumentos!")
  if(class(spec) != "uGARCHspec") stop("roll_fit: spec deve ser da classe uGARCspec")
  
  tmp.list <- map(models, 
                  ~switch (.x,
                           cevt = roll_fit_cevt(data, window.size, n.roll, spec),
                           uevt = roll_fit_uevt(data, window.size, n.roll, spec),
                           unorm = roll_fit_unorm(data, window.size, n.roll),
                           ut = roll_fit_ut(data, window.size, n.roll),
                           riskmetrics = roll_fit_riskmetrics(data, window.size, n.roll),
                           cat("Modelo não definido:", .x)
                           ) # fim switch
                  ) # fim map
  names(tmp.list) <- models
  tmp.list <- enframe(tmp.list)
  colnames(tmp.list) <- c("model_type", "roll")
  tmp.list <- subset(tmp.list, subset = !is.null(tmp.list$roll))
  cat("\nroll_fit do ativo fim:", as.character(Sys.time()))
  return(tmp.list)
}
# roll_fit_cevt -----------------------------------------------------------
roll_fit_cevt <- function(data, window.size, n.roll, spec) {
  tic <- Sys.time()
  # Check para os argumentos
  if(!is.xts(data)) stop("roll_fit_cevt: data deve ser um xts")
  if(any(is.null(window.size), is.null(n.roll), is.null(spec)))
    stop("roll_fit_cevt: Devem ser passados todos os argumentos!")
  if(class(spec) != "uGARCHspec") stop("roll_fit_cevt: spec deve ser da classe uGARCspec")
  
  # Deve-se fazer o fit
  # Extrair residuos padronizados
  # Aplicar o gpfFit
  # Calcular zq e sq
  # Para cada interacao. n.roll vezes
  # Devolver tudo em um data.frame
  # data deve ser o xts com as perdas de todo o periodo
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
      # Zq = mu_t+sqrt(sigma_t)*zq
      # Sq = mu_t+sqrt(sigma_t)*sq
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
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("beta1", "xi", "beta", "xi_se", "beta_se")]
  toc <- Sys.time()
  cat("\ncevt:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
  # Os valores retornados das medidas de risco devem ser comparadas
  # com os valores realizados NO DIA SEGUINTE a data onde foram calculadas
} # fim da roll_fit_cevt


# roll_fit_unorm ----------------------------------------------------------
# Ajusta os dados para um modelo Normal incondicional
roll_fit_unorm <- function(data, window.size, n.roll) {
  tic <- Sys.time()
  # Check para os argumentos
  if(!is.xts(data)) stop("roll_fit_unorm: data deve ser um xts")
  if(any(is.null(window.size), is.null(n.roll)))
    stop("roll_fit_unorm: Devem ser passados todos os argumentos!")
  
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
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("mu", "sigma")]
  toc <- Sys.time()
  cat("\nunorm:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
}

# roll_fit_ut ----------------------------------------------------------
# Ajusta os dados para um modelo t-Student incondicional
roll_fit_ut <- function(data, window.size, n.roll) {
  tic <- Sys.time()
  # Check para os argumentos
  if(!is.xts(data)) stop("roll_fit_ut: data deve ser um xts")
  if(any(is.null(window.size), is.null(n.roll)))
    stop("roll_fit_ut: Devem ser passados todos os argumentos!")
  
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
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("mu", "sigma", "shape")]
  toc <- Sys.time()
  cat("\nut:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
}

# roll_fit_uevt -----------------------------------------------------------
roll_fit_uevt <- function(data, window.size, n.roll, spec) {
  tic <- Sys.time()
  # Check para os argumentos
  if(!is.xts(data)) stop("roll_fit_uevt: data deve ser um xts")
  if(any(is.null(window.size), is.null(n.roll), is.null(spec)))
    stop("roll_fit_uevt: Devem ser passados todos os argumentos!")
  if(class(spec) != "uGARCHspec") stop("roll_fit: spec deve ser da classe uGARCspec")
  
  # Deve-se fazer o fit
  # Extrair residuos padronizados
  # Aplicar o gpfFit
  # Calcular zq e sq com base na media e variancia INCONDICIONAIS
  # Para cada interacao. n.roll vezes
  # Devolver tudo em um data.frame (ou lista)
  # data deve ser o xts com as perdas de todo o periodo
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
      mu_t <- uncmean(garch.fit)
      sigma_t <- sqrt(uncvariance(garch.fit)) 
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
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("beta1", "xi", "beta", "xi_se", "beta_se")]
  toc <- Sys.time()
  cat("\nuevt:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
  # Os valores retornados das medidas de risco devem ser comparadas
  # com os valores realizados NO DIA SEGUINTE a data onde foram calculadas
} # fim da roll_fit_uevt

# roll_fit_riskmetrics -----------------------------------------------------------
roll_fit_riskmetrics <- function(data, window.size, n.roll){
  tic <- Sys.time()
  # Check para os argumentos
  if(!is.xts(data)) stop("roll_fit_riskmetrics: data deve ser um xts")
  if(any(is.null(window.size), is.null(n.roll)))
    stop("roll_fit_riskmetrics: Devem ser passados todos os argumentos!")
  # Eh uma especificacao de Garch(1,1) com mu = 0, omega = 0, lambda = beta1 e 1-lambda = alpha1
  # Como os parametros do Garch sao fixos, pode-se utilizar o metodo ugarchfilter
  lambda <- 0.94 # Valor apontado como ideal para dados diarios
  ruspec <- ugarchspec(mean.model = list(armaOrder = c(0,0),
                                         include.mean = FALSE),
                       variance.model = list(model = "iGARCH",
                                             garchOrder = c(1, 1)),
                       distribution.model = "norm",
                       fixed.pars = list(omega = 0,
                                         alpha1 = (1-lambda))) # Beta eh calculado no modelo iGarch
  filter <- ugarchfilter(ruspec, data[(window.size+1):(window.size+n.roll)])
  resid <- residuals(filter)
  sigma <- sigma(filter) # variancia!! eh necessario tirar a raiz para obter o desv. padrao
  Zq975 <- sigma*qnorm(0.975)
  Zq990 <- sigma*qnorm(0.990)
  Sq975 <- (coredata(sigma)*dnorm(qnorm(0.975)))/0.025 # Eq 4.7 p. 37 de Pfaff2013
  Sq990 <- (coredata(sigma)*dnorm(qnorm(0.990)))/0.01
  
  ans <- cbind(lambda, Zq975, Zq990, Sq975, Sq990) # 
  ans <- ans[-dim(ans)[1],]           
  ans <- xts(ans, order.by = index(data[(window.size+2):(window.size+n.roll)]))
  colnames(ans) <- c("lambda", "Zq975", "Zq990", "Sq975", "Sq990") # 
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  risk <- tibble(coverage = c(0.025, 0.01),
                 VaR.xts = list(ans$Zq975, ans$Zq990),
                 ES.xts = list(ans$Sq975, ans$Sq990))
  param <- ans[, c("lambda")]
  toc <- Sys.time()
  cat("\nriskmetrics:", toc-tic, attr(toc-tic, which = "units"))
  return(list(risk.tbl = risk, param.xts = param))
}

# es_test -----------------------------------------------------------------
# Teste nao parametrico para os residuos das violacoes ao VaR, conforme 
# teste de ES de MacNeil2000
# ATENCAO!! Este teste so faz sentido ser realizados apos os teste de VaR
# aprovarem o modelo
es_test <- function(alpha = 0.0275, losses, ES, VaR, conf.level = 0.95, n.boot = 1000) {
  # Check for univariate series
  if(!all(dim(losses)[2] == 1,
          dim(ES)[2] == 1,
          dim(VaR)[2] == 1))
    stop("\nSeries are not univariate!")
  # Check for the same length in all series
  N <-  length(losses)
  if(!all(length(VaR) == N,
          length(ES) == N))
    stop("\nLength of series are not equal!")
  idx <-  which(losses > VaR)
  s <-  ES[idx] # ES values when VaR violation occurred
  x <- losses[idx] # losses that violated VaR
  
  # One-sided test. H0: mean of x is less than or equal to mean of s
  # We do not want to reject H0
  boot <- boot.t.test(x, s, reps = n.boot, alternative = "greater")
  
  ans <-  tibble(expected.exceed = floor(alpha*N),
                 actual.exceed = length(idx)) %>% 
    bind_cols(boot)
  # conditional expected shortfall is systematically underestimated
  ans$H0 <- "Mean of Excess Violations of VaR is less than or equal to zero"
  ans$Decision <- ifelse(ans$p.value < (1 - conf.level), "Reject H0", "Failed to reject H0")
  return(ans)
} # end of es_test

# boot.t.test -------------------------------------------------------------
# Codigo importado de https://github.com/tpepler/nonpar
# Adaptado apenas para as checagens dos argumentos
boot.t.test <- function(x, y, reps = 1000, mu = 0, alternative = c("two.sided", "less", "greater")){
  # Bootstrap t-test as described in Efron and Tibshirani (1993), (Algorithm 16.2, p224)
  if(is.null(x) | is.null(y)) 
    stop("\nArguments to boot.t.test cannot be NULL!")
  if((length(x) <= 1) | (length(y) <= 1)) 
    stop("\nboot.t.test: Sample lengths cannot be less than 10!")
  
  nx <- length(x)
  ny <- length(y)
  t.obs <- (mean(x) - mean(y) - mu) / sqrt(var(x) / nx + var(y) / ny)
  comb.mean <- mean(c(x, y))
  x.c <- x - mean(x) + comb.mean
  y.c <- y - mean(y) + comb.mean
  t.boot <- rep(NA, times = reps)
  
  bootFunc <- function(){
    bootx <- x.c[sample(1:nx, size = nx, replace = TRUE)]
    booty <- y.c[sample(1:ny, size = ny, replace = TRUE)]
    return((mean(bootx) - mean(booty) - mu) / sqrt(var(bootx) / nx + var(booty) / ny))
  }
  
  t.boot <- replicate(reps, expr = bootFunc())
  
  if(alternative[1] == "two.sided"){
    pval <- length(t.boot[abs(t.boot) >= abs(t.obs)]) / reps
    h1phrase <- "not equal to"
  }
  
  if(alternative[1] == "less"){
    pval <- length(t.boot[t.boot <= t.obs]) / reps
    h1phrase <- "less than"
  }
  
  if(alternative[1] == "greater"){
    pval <- length(t.boot[t.boot >= t.obs]) / reps
    h1phrase <- "greater than"
  }
  
  cat("\nBootstrap Two Sample t-test\n")
  cat(paste("\nt = ", round(t.obs, 3), ", p-value = ", round(pval, 4), "\n", sep=""))
  cat(paste("Alternative hypothesis: true difference in means is ", h1phrase, " ", mu, "\n\n", sep=""))
  
  return(tibble(mu0 = mu,
                statistic = t.obs,
                alternative = alternative[1],
                p.value = pval))
}

