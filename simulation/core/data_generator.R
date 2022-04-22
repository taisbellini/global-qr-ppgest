# Simulating the data

data_gen = function(N = 100, D=30, D_signal=20, D_strong=9){
  
  # X = matrix(rnorm(N*D, mean = 10),N,D)
  X = matrix(runif(N*D),N,D)
  X[,1] = 1
  
  U = runif(N)
  
  Y = rep(0,N)
  for (n in 1:N){
    Y[n] = t(X[n,])%*%beta(U[n], D, D_signal, D_strong) # generating Y via the Fundamental Theorem of Simulation
  }
  
  return (list(
    "X" = X,
    "y" = Y,
    "U" = U
  ))
}

beta = function(tau, D, D_signal, D_strong){
  D_weak = D_signal - D_strong
  Dnoise = D - D_signal
  
  coef_strong = 2/D_signal
  coef_weak = .1/D_signal
  theta = c(rep(coef_strong,D_strong),rep(coef_weak,D_weak),rep(0,Dnoise))
  sapply(1:D, function(d) theta[d]*tau^(d - 1))
}

data_gen_strong = function(N = 100, D=20, D_signal=10, D_strong=7){
  
  # X = matrix(rnorm(N*D, mean = 10),N,D)
  X = matrix(runif(N*D),N,D)
  X[,1] = 1
  
  U = runif(N)
  
  Y = rep(0,N)
  for (n in 1:N){
    Y[n] = t(X[n,])%*%beta_strong(U[n], D, D_signal, D_strong) # generating Y via the Fundamental Theorem of Simulation
  }
  
  return (list(
    "X" = X,
    "y" = Y,
    "U" = U
  ))
}

beta_strong = function(tau, D, D_signal, D_strong){
  D_weak = D_signal - D_strong
  Dnoise = D - D_signal
  
  coef_strong = 2/1
  coef_weak = .1/1
  theta = c(rep(coef_strong,D_strong),rep(coef_weak,D_weak),rep(0,Dnoise))
  sapply(1:D, function(d) theta[d]*tau^(d - 1))
}




