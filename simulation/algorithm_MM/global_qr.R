#### GLOBAL QR ####

## Objective function to minimize ##
R_vec = function(a_vec, w, zeta, y, X, tau.grid, phi, eta, f = "lasso", lambda = 0){
  L = nrow(phi)
  N = length(y)
  M = length(tau.grid)
  R.eval = 0
  
  if (f == "lasso"){
    penalty = (lambda/2)*sum((c(w)*((a_vec^2/c(zeta))+c(zeta))))
  } else {
    w = matrix(rep(w, L), nrow= L, byrow = T)
    zeta = matrix(rep(zeta, L), nrow= L, byrow = T)
    a = matrix(a_vec^2, ncol=L, byrow=FALSE)
    a_sum = apply(a, 1, sum)
    a_glasso = sqrt(a_sum)
    penalty = (lambda/2)*sum((c(w)*((a_glasso^2/c(zeta))+c(zeta))))
  }
  
  for (n in 1:N){
    for (m in 1:M){
      res = t(c(X[n,]%*%t(phi[,m])))%*%a_vec
      res = res[1,1]
      R = (2*tau.grid[m] - 1)*(y[n] - res) + 
        (1/2)*( ((y[n] - res)^2 /eta[n,m])+ eta[n,m])
      R.eval = R.eval + R
    }
  }
  return(R.eval + penalty)
}

## Calculate a_vec ##
## TODO: for's just once
a_vec_opt = function(y, w, X, phi, tau.grid, eta, zeta, f = "lasso", lambda = 0){
  N = nrow(X)
  M = length(tau.grid)
  L = nrow(phi)
  
  # Define penalty term
  if (f == "lasso"){
    penalty = lambda*diag(c(w/zeta))
  } else {
    w = matrix(rep(w, L), nrow= L, byrow = T)
    zeta = matrix(rep(zeta, L), nrow= L, byrow = T)
    penalty = lambda*diag(c(w/zeta))
  }

  el1 = 0
  el2 = 0
  for (n in 1:N){
    for (m in 1:M){
      Xn_phim = Xn_phim_vec(X[n,], phi[,m])
      el1 = el1 - ((1/eta[n,m])*(Xn_phim%*%t(Xn_phim)))
    }
  }
  el1 = el1 - penalty
  
  #solve_el1 = tryCatch({
  #  el1 = solve(el1)
  #}, warning = function(w){
  #  print(paste("Error: ", w))
  #  el1 = NULL
  #}, error = function(e){
  #  print(paste("Error: ", e))
  #  el1 = NULL
  #}, finally = function(f){
  #  el1 = el1
  #})
  
  el1 = solve(el1)
  #if(is.null(el1)) return("error")
  
  for (m in 1:M){
    for (n in 1:N){
      Xn_phim = Xn_phim_vec(X[n,], phi[,m])
      el2 = el2 + ((1-(2*tau.grid[m]))-(y[n]/eta[n,m]))*Xn_phim
    }
  }
  return (el1%*%el2)
}

a_vec_opt2 = function(y, w, X, phi, tau.grid, eta, zeta, f = "lasso", lambda = 0){
  N = nrow(X)
  M = length(tau.grid)
  L = nrow(phi)
  
  # Define penalty term
  if (f == "lasso"){
    penalty = lambda*diag(c(w/zeta))
  } else {
    w = matrix(rep(w, L), nrow= L, byrow = T)
    zeta = matrix(rep(zeta, L), nrow= L, byrow = T)
    penalty = lambda*diag(c(w/zeta))
  }
  
  el1 = 0
  el2 = 0
  for (n in 1:N){
    for (m in 1:M){
      Xn_phim = Xn_phim_vec(X[n,], phi[,m])
      el1 = el1 - ((1/eta[n,m])*(Xn_phim%*%t(Xn_phim)))
    }
  }
  el1 = diag(sqrt(c(zeta)))%*%el1%*%diag(sqrt(c(zeta)))
  el1 = el1 + lambda*diag(rep(1, D*L))
  el1 = solve(el1)
  el1 = diag(sqrt(c(zeta)))%*%el1%*%diag(sqrt(c(zeta)))
  for (m in 1:M){
    for (n in 1:N){
      Xn_phim = Xn_phim_vec(X[n,], phi[,m])
      el2 = el2 + ((1-(2*tau.grid[m]))-(y[n]/eta[n,m]))*Xn_phim
    }
  }
  return (el1%*%el2)
}

## Calculate eta and zeta ## 

eta_gen = function(y, X, a_vec, phi, pi_p = 0.01){
  N = length(y)
  M = ncol(phi)
  eta_matrix = matrix(rep(0, N*M), ncol = M)
  a = matrix(a_vec, ncol=L, byrow=FALSE)
  for (n in 1:N){
    for (m in 1:M){
      eta_matrix[n,m] = sqrt(abs(y[n] - X[n,]%*%a%*%phi[,m])^2 + pi_p)
    }
  }
  return(eta_matrix)
}

zeta_gen = function(a, f = "lasso", pi_p = 0.01){
  
  if (f == "lasso"){
    zeta = sqrt(a^2 + pi_p)
  }
  else if (f == "gLasso"){
    a = a^2
    a_sum = apply(a, 1, sum)
    zeta = sqrt(a_sum)
    # Precision parameter adjustment
    zeta_p = sapply(zeta, function(value) {max(value, pi_p)})
  } else {
    zeta = matrix(1, nrow(a), ncol(a))
  }
  
  return(zeta)
}


## Aux functions ##

Xn_phim_vec = function(Xn, phim){
  return(c(Xn%*%t(phim)))
}

check_function = function(u, tau){
  rho = (1/2)*((2*tau - 1)*u + abs(u))
  return(rho)
}

R_rho = function(y, X, a_k, phi, tau.grid){
  
  N = nrow(X)
  M = length(tau.grid)
  
  R = 0
  for(n in 1:N){
    for (m in 1:M){
      R = R + check_function(y[n] - t(X[n,])%*%a_k%*%phi[,m], tau.grid[m])
    }
  }
  return (R)
}

lambda_opt = function(a_0, phi, X, y, we, tau.grid = seq(from = .1, to = .9, by = .1), f = "lasso", lambdas_to_try = 10^seq(-3, 2, length.out = 20)){
  
  bic_vec = c()
  aic_vec = c()
  
  for (lambda_ind in seq(lambdas_to_try)){
    
    mod = global_qr_opt(a_0, y, we, X, phi, tau.grid, f = f, lambda = lambdas_to_try[lambda_ind])
    
    if(typeof(mod$bhat) == "character" & mod$bhat == "error") {
      bic_vec[lambda_ind] = NA
      aic_vec[lambda_ind] = NA
      next
    }
    
    mod_coefs = mod$bhat
    df = nrow(mod_coefs[rowSums(mod_coefs) != 0,])
    if (is.null(df)) df = 1
    #dev = mod$R_k
    dev = R_rho(y, X, mod$ahat, phi, tau.grid)
    n = nrow(X)
    
    aic = log(dev) + (n^(-1))*df
    bic = log(dev) + ((2*n)^(-1))*log(n)*df*1
    
    aic_vec[lambda_ind] = aic
    bic_vec[lambda_ind] = bic
    
  }
  print(aic_vec)
  print(bic_vec)
  return(
    list(
      "aic" = lambdas_to_try[which.min(aic_vec)],
      "bic" = lambdas_to_try[which.min(bic_vec)]
    ))
}

convergence = function(R_old, R_k, tol=1e-6){
  if (abs((R_old - R_k)/R_old) < tol) {return(T)}
  return(F)
}

convergence2 = function(a_old, a_k, tol=1e-6){
  if (max(abs((a_old - a_k))) < tol) {return(T)}
  return(F)
}

## Optimization procedure - core ##

global_qr_opt = function(a_0, y, we, X, phi, tau.grid, f = "lasso", lambda = 0){
  
  # Get parameters
  M = length(tau.grid)
  L = nrow(phi)
  D = ncol(X)
  N = nrow(X)
  
  # Tolerance to consider coefficient equal zero
  tol = 1e-6
  
  # Initialize algorithm
  a_vec_0 = c(a_0) # vec(a)
  eta_0 = eta_gen(y, X, a_vec_0, phi) #eta matrix
  zeta_0 = zeta_gen(a_0, f = f) # zeta matrix if lasso, zeta vec if gLasso
  
  a_vec_k = a_vec_0
  eta_k = eta_0
  zeta_k = zeta_0
  
  R_old = R_vec(a_vec_0, we, zeta_0, y, X, tau.grid, phi, eta_0, f = f, lambda = lambda)
  a_old = a_0
  
  maxit = 1000
  
  for(i in 1:maxit){
    a_vec_k = a_vec_opt(y, we, X, phi, tau.grid, eta_k, zeta_k, f = f, lambda = lambda)
    if(typeof(a_vec_k) == "character") return ("error")
    eta_k = eta_gen(y, X, a_vec_k, phi)
    a_k = matrix(a_vec_k, nrow = D)
    zeta_k = zeta_gen(a_k, f = f)
    
    R_k = R_vec(a_vec_k, we, zeta_k, y, X, tau.grid, phi, eta_k, f = f, lambda = lambda)
    
    if(R_k > R_old){
      print(paste("subiu: ", R_k-R_old, "iteração: ", i))
      print(lambda)
      
    }
    
    if (convergence2(a_old, a_k)) break
    
    R_old = R_k
    a_old = a_k
    
  }
  
  ahat = matrix(a_vec_k, D, L)
  ahat_tol = apply(ahat, 1, function(row) {
    if (sum(abs(row))< tol){
      return(rep(0, length(row)))
    }
    return(row)
  })
  bhat = t(ahat_tol)%*%phi
  
  return(list(
    "ahat" = ahat,
    "bhat" = bhat,
    "ahat_tol" = t(ahat_tol), 
    "R_k" = R_k
  ))
}


## Optimization procedure - root ## 

global_qr = function(a_0, tau.grid = c(0.5), phi = matrix(), X = matrix(), y = c(), lambda = 0, f = "lasso", w = F, criteria = "bic"){
  
  # Calculate weights
  if(w){
    we = weights(a, Y, X, tau.grid, phi)
  } else {
    if (f == "lasso"){
      we = matrix(1, D, L)
    } else {
      we = rep(1, D)
    }
  }
  
  if(typeof(we) == "character") return (list("bhat" = "error"))
  
  if(lambda == "opt"){
    opt_lambda = lambda_opt(a_0, phi, X, y, we, tau.grid = tau.grid, f = f)
    if (criteria == "bic") {
      lambda = opt_lambda$bic
    } else {
      lambda = opt_lambda$aic
    }
  }
  
  opt_global_qr = global_qr_opt(a_0, y, we, X, phi, tau.grid, f = f, lambda = lambda)
  
  if(typeof(we) == "character") return (list("bhat" = "error"))
  
  return(list(
    "ahat" = opt_global_qr$ahat,
    "bhat" = opt_global_qr$bhat,
    "ahat_tol" = opt_global_qr$ahat_tol,
    "lambda" = lambda
  ))
}




