#### CVXR 

## aux functions #####
lambda_opt_CVXR = function(phi, X, y, we, tau.grid = seq(from = .1, to = .9, by = .1), f = "lasso", lambdas_to_try = 10^seq(-3, 3, length.out = 100), s=2, type = "par"){
  
  tol = 1e-6
  
  D = ncol(X)
  L = nrow(phi)
  
  ahat_vec = array(0, c( D, L, length(lambdas_to_try)))
  
  for(lambda_ind in seq(lambdas_to_try)){
  
    ahat = global_qr_CVXR_opt(we, taus = tau.grid, phi = phi, X = X, y = y, lambda = lambdas_to_try[lambda_ind], f = f, s=s, type= type)
    
    if(typeof(ahat) == "character"){
      ahat_vec[,,lambda_ind] = rep(0, D*L)
    } else {
      
      ahat_tol = ahat
      ahat_tol[abs(ahat_tol)<tol] = 0
      rs = rowSums(ahat_tol)
      
      ahat_vec[,,lambda_ind] = ahat_tol
      
      if ( length(rs[rs>0]) <= 1 ){
        return(ahat_vec)
      }
    }
  }
  
  return(ahat_vec)
}

R_CVXR = function(a, Y, X, TAUS, phi, type = "par"){
  if(type == "naive"){
    B = a
  }
  else {
    B = a%*%phi # candidate optimizer
  }
  R.eval = (2*TAUS-1)*(Y-X%*%B) + abs(Y-X%*%B)
  return(sum(R.eval)) # mean of the preceding objective functions
}

### Save ahat and bhat found

penalty_CVXR = function(a, f, w, s = 2){
  
  if(f == "lasso"){
    penalty = sum(w*abs(a))
  } else {
    #a_norm = apply(value(a), 1, function(row){(sum(abs(row)^s))^(1/s)})
    penalty = w[1]*p_norm(a[1,], s)
    for (d in 2:nrow(a)){
      penalty  = penalty + w[d]*p_norm(a[d,], s)
    }
  }
  return(penalty)
}

# sum(sapply(1:D, function(d) sqrt(sum(a[d,]^2))))
weights_CVXR = function(y, X, tau.grid, phi, f, s=2, type = "par", w = F){
  
  we = numeric()
  p = 1
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


## Optimization procedure - core ####
global_qr_CVXR_opt = function(we, taus = c(0.5), phi = matrix(), X = matrix(), y = c(), lambda = 0, f = "lasso", s = 2, type = "par"){
  
  # Get parameters
  M = length(taus)
  L = nrow(phi)
  D = ncol(X)
  N = nrow(X)
  Y = matrix(rep(y,M),N,M)
  TAUS = matrix(rep(taus,N),N,M, byrow = TRUE)
  

  a = Variable(D,L)
  
  objective = R_CVXR(a, Y, X, TAUS, phi, type) + lambda*penalty_CVXR(a, f, we, s = s)
  problem = Problem(Minimize(objective))
  result = solve(problem)
  
  if(result$status == 'solver_error') {
    return ("error")
  }else if(is.na(result$getValue(a))) {
    return("error")
  }else {
    ahat = result$getValue(a)
    return(ahat)
  }
}

## Optimization procedure ####
global_qr_CVXR = function(taus = c(0.5), phi = matrix(), X = matrix(), y = c(), lambda = 0, f = "lasso", s = 2, w = F, lambdas = 10^seq(-3, 3, length.out = 100), type = "par"){
  
  D = ncol(X)
  L = nrow(phi)

  we = weights_CVXR(y, X, taus, phi, f, s = s, type = type, w = w)
  
  
  if (lambda == "opt"){
    opt_global_qr_CVXR = lambda_opt_CVXR(phi, X, y, we, tau.grid = taus, f = f, lambdas_to_try = lambdas, s= s, type = type)
  } else{
    opt_global_qr_CVXR = global_qr_CVXR_opt(we, taus = taus, phi = phi, X = X, y = y, lambda = lambda, f = f, s= s, type = type)
    
    return(opt_global_qr_CVXR)
  }
}


