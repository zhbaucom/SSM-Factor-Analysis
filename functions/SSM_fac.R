SSM_fac <- function(
  data, outcome_vars, id_var, time_var, G_form,
  its = 10000, burn = floor(its/2), wait = floor(its/4),
  mu_G = matrix(0, length(outcome_vars), ncol(G_form)), var_G = 10, #mean and variance priors for G
  c0 = 0.01, d0 =0.01, #priors for sigeps2
  m0 = 0, C0 = 10, #priors for alpha
  seed = NULL
){
  
  if(!is.null(seed))set.seed(seed)
  K <- length(outcome_vars)
  Q <- ncol(G_form)
  
  
  #Initialize y and time lists
  y <- list()
  time <- list()
  alpha.star <- list()
  #Create list to hold each y and time
  uid <- unique(data[[id_var]])
  nid <- N <- length(uid)
  Neta <- nrow(data) - nid
  for(i in 1:nid){
    tfid <- data[id_var] == uid[i]
    y[[i]] <- t(as.matrix(FULL[tfid, outcome_vars]))
    time[[i]] <- data[[time_var]][tfid]
    alpha.star[[i]] <- matrix(m0, nrow = nrow(y[[i]]), ncol = ncol(y[[i]]))
  }
  
  time_diff <- map(time, diff)
  
  # Initialize Variances
  V <- diag(1, K, K)
  W <- diag(1, Q, Q)
  

  
  #Variance Holders  
  vcovWish <- matrix(NA, its, K)
  wcovWish <- array(NA, dim = c(Q, Q, its))
  
  ###Priors for ETA
  prior.nu.eta <- K - 1 + 0.01
  ## prior.nu.eta <- cs + 1
  prior.Gamma.eta <- diag(0.01, K)
  ## prior.Gamma.eta <- diag(1, cs)
  
  

  

  
  G.star <- mu_G
  # G.star[1,1] <- 0
  
  G.track <- array(dim = c(dim(G_form), its))
  

  
  ybig <- t(as.matrix(data[outcome_vars]))
  alpha.track <- list()
  
  for(i in 1:its){
    
    ### G transofrmed Alphas
    alpha.star <- alpha.track[[i]] <- lapply(1:nid, function(yi){
      ffbs.joint(
        y = y[[yi]], G = G.star, 
        V = V, W = W, 
        m0 = 0, C0 = diag(10, Q, Q), 
        timeDiff = time_diff[[yi]]
      )$x
    })
    
    
    ### Posterior of G
    abig <- do.call("cbind", alpha.star)
    
    
    aybig <- tcrossprod(ybig, abig)
    
    ata.inv <-  solve(tcrossprod(abig)/sigeps2 + diag(1/var_G, Q, Q))
    mu.g <- (mu_G/var_G + aybig/sigeps2) %*% ata.inv
    
    
    G.track[,,i] <- G.star <- abs((lapply(1:K, function(k){
      mvtnorm::rmvnorm(1, mean = mu.g[k,], sigma = ata.inv)
    }) %>%
      do.call("rbind", .)) * G_form)
    
    
    
    ### Variances
    if(i > wait){
      fp <- lapply(1:N, function(i){
        apply(alpha.star[[i]], 1, function(z)diff(z)) / sqrt( time_diff[[i]] )
        
      }) %>%
        do.call("rbind", .) %>%
        crossprod()
      
      
      W <- cov2cor(cIRT::riwishart(Neta, fp))
      
      wcovWish[,,i] <- W
    }
    

    
    ###### Sigma eps
    
    or <- lapply(1:N, function(i){
      (y[[i]] - G.star %*% alpha.star[[i]])^2
    })  %>%
      do.call("cbind", .) %>%
      apply(1, function(x)sum((x)))
    
    
    sigma2.eps.star <- sapply(or, function(x)1/rgamma(1, ((N*TT)/2 +c0), d0+x/2))
    vcovWish[i,] <- sigma2.eps.star
    V <- sigma2.eps.star * diag(K)
    
    if((i %% 10) == 0)replaceMessage(i)
    
  }
  
  list(G = G.track, alpha = alpha.track, Eps = vcovWish, Eta = wcovWish)
  
}