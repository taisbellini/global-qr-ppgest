#### Utils dissertacao ####

#### Helper functions ####

# Generate the phi matrix based on the taus grid and number of basis vectors
phi_generator = function(L, taus){
  M = length(taus)
  phi = matrix(0,L,M) # matrix to store basis vectors
  colnames(phi) = as.character(taus)
  
  # finite dimensional analogue of Legendre polynomials
  phi[1,] = rep(1,M)
  phi[1,] = phi[1,]/sqrt(sum(phi[1,]^2))
  phi[2,] = lm(taus~0+phi[1,])$res
  phi[2,] = phi[2,]/sqrt(sum(phi[2,]^2))
  if (L>=3) {
    for (ell in 3:L){
      phi[ell,] = lm(taus^(ell-1)~t(phi[2:ell,]))$res
      phi[ell,] = phi[ell,]/sqrt(sum(phi[ell,]^2))
    }
  }
  return(phi)
}

tau.grid_generator = function(L, M, delta = .1){
  Lgen = L-1
  Mgen = M-1 # damn you, R! should allow indexing matrices from 0
  f = function(x) .5 + (1-2*delta)*x/2
  tau.grid = f(cos(pi*0:Lgen/Lgen))
}

phi_generator2 = function(L, tau.grid, delta = .1){
  
  f.inv = function(t) 1/(2*delta-1) + 2*t/(1-2*delta)
  
  phi.fun = function(tau, m){
    acos_r = ifelse(tau == delta, acos(-1), ifelse(tau == (1-delta), acos(1), acos(f.inv(tau))))
    return(cos(m*acos_r))
  }
  
  phi.matrix = matrix(0,L,L)
  
  for (ell in 0:(L-1)){
    phi.matrix[ell+1,] = phi.fun(tau.grid,ell)
  }
  
  L = nrow(phi.matrix) # "true" L
  M = ncol(phi.matrix) # "true" M
  
  return (phi.matrix)
}

rho_function = function(u, t){
  rho = (1/2)*((2*t-1)*u + abs(u))
  return(rho)
}


#### Aux functions to find AIC/BIC ####

R_function = function(ahat, Y, X, TAUS, phi, type = "par"){
  if(type == "naive"){
    B = ahat
  }else {
    B = ahat%*%phi # candidate optimizer
  }
  R.eval = (2*TAUS-1)*(Y-(X%*%B)) + abs(Y-(X%*%B))
  return(sum(R.eval)) # mean of the preceding objective functions
}

weights_function = function(y, X, tau.grid, phi, f, w, type){
  we = numeric()
  p = 1
  s = 2
  D = ncol(X)
  L = nrow(phi)
  
  if(w){ 
    a_rq = rq_opt(y, X, tau.grid, phi, type=type)
    
    if (f == "lasso"){ #adaLasso
      we = (abs(a_rq))^(-p)
    }else { #group adaLasso
      a_norm = apply(a_rq, 1, function(row){(sum(abs(row)^s))^(1/s)})
      we = (a_norm)^(-p)
    }
  } else { 
    if(f == "lasso"){ #lasso
      we = matrix(1, nrow = D, ncol = L)  
    } else if (f == "gLasso"){ #group Lasso
      we = rep(1, D)
    }
  }
  
  return(we)
}

penalty_function = function(a, f, w, s = 2){
  
  if(f == "lasso"){
    penalty = sum(w*abs(a))
  } else {
    a_norm = apply(a, 1, function(row){(sum(abs(row)^s))^(1/s)})
    penalty = sum(w*a_norm)
  }
  return(penalty)
}

## canonical estimator ####

rq_opt = function(y, X, tau.grid, phi, type = "par"){
  X = X[,2:ncol(X)]
  qrfit = rq(y ~ X, tau = tau.grid)
  B_QR = coef(qrfit)
  
  if (type == "naive"){
    return(B_QR)
  }else{
    return (B_QR %*% solve(phi))
  }
}

