source("../core/data_generator.R")
source("../core/utils.R")
source("../core/global_qr_CVXR.R")
source("global_qr.R")


convergence = function(R_old, R_k, tol=1e-15){
  if (abs((R_old - R_k)/R_old) < tol) {return(T)}
  return(F)
}

## vec(a) in optim function
# sugestão
R_vec = function(b, y, X, tau.grid, phi, eta, d_index, l_index, pi=.001){
  L = nrow(phi)
  N = length(y)
  M = length(tau.grid)
  R.eval = 0
  a = matrix(b, ncol=L, byrow=FALSE)
  for (n in 1:N){
    for (m in 1:M){
      res = (t(X[n,])%*%a%*%phi[,m])[1,1]
      R = (2*tau.grid[m] - 1)*(y[n] - res) + 
        (1/2)*( ((y[n] - res)^2 /max(pi,eta[n,m]))+ eta[n,m] )
      R.eval = R.eval + R
    }
  }
  return(R.eval)
}

eta_vec_gen = function(y, X, a_vec, phi){
  N = length(y)
  M = ncol(phi)
  eta_matrix = matrix(rep(0, N*M), ncol = M)
  a = matrix(a_vec, ncol=L, byrow=FALSE)
  for (n in 1:N){
    for (m in 1:M){
      eta_matrix[n,m] = abs(y[n] - t(X[n,])%*%a%*%phi[,m])
    }
  }
  return(eta_matrix)
}

## todo: for's just once
a_vec_opt = function(y, X, phi, tau.grid, eta){
  N = nrow(X)
  M = length(tau.grid)
  
  el1 = 0
  el2 = 0
  for (n in 1:N){
    for (m in 1:M){
      Xn_phim = Xn_phim_vec(X[n,], phi[,m])
      el1 = el1 - (1/eta[n,m])*(Xn_phim%*%t(Xn_phim))
    }
  }
  for (m in 1:M){
    for (n in 1:N){
      Xn_phim = Xn_phim_vec(X[n,], phi[,m])
      el2 = el2 + ((1-(2*tau.grid[m]))-(y[n]/eta[n,m]))*Xn_phim
    }
  }
  return (solve(el1)%*%el2)
}

Xn_phim_vec = function(Xn, phim){
  return(as.vector(Xn%*%t(phim)))
}

## Gradient function

gr = function(v, y, X, tau.grid, phi, eta, d_index, l_index, pi = 0){
  gradient_v = vector()
  M = length(tau.grid)
  N = nrow(X)
  D = ncol(X)
  L = nrow(phi)
  J = length(v)
  
  a = matrix(v, ncol=L, byrow=FALSE)
  
  for (j in 1:J){
    Rj = 0
    for(n in 1:N){
      for(m in 1:M){
        Rj = Rj + (1-2*tau.grid[m])*(X[n,d_index[j]]*phi[l_index[j],m]) - ((y[n] - t(X[n,])%*%a%*%phi[,m]) / max(pi, eta[n,m]))*(X[n,d_index[j]]*phi[l_index[j],m])
      }
    }
    gradient_v[j] = Rj
  }
  return(gradient_v)
}


gr2 = function(v, y, X, tau.grid, phi, eta, pi = 0){
  M = length(tau.grid)
  N = nrow(X)
  D = ncol(X)
  L = nrow(phi)
  
  R = 0
  for(n in 1:N){
    for(m in 1:M){
      Xn_phim = Xn_phim_vec(X[n,], phi[,m])
      R = R + ((1-(2*tau.grid[m]))*Xn_phim - ((y[n] - t(Xn_phim)%*%v) / max(pi, eta[n,m])))*Xn_phim
    }
  }
  return(R)
}

gr_mat = function(a, Y, X, phi, tau, eta){
  foo = matrix(0, D, L)
  for (n in 1:N){
    for (m in 1:M){
      Xnloc = matrix(X[n,], D, 1)
      phimloc = matrix(phi[,m], L, 1)
      foo = foo + c((1-2*tau[m]) - (Y[n] - t(Xnloc)%*%a%*%phimloc)/eta[n,m])*(Xnloc%*%t(phimloc))
    }
  }
  return(foo)
}


#### data generation

source("data_generator.R")
source("utils.R")
### Simulate cross section Data ####
n = 1000
N = n
epsilon = rnorm(n, sd=.01)
Dgen = 5 # total vars
Dsignal = 3 # vars = 0

# data gen function from data_generator.R
data = data_gen(epsilon, n = n, D=Dgen, Dsignal = Dsignal)
X = data$X
X_1 = cbind(1, X)
X = X_1
y = data$y[,1]

# params def
tau.grid = seq(from = .1, to = .9, by = .1)
#tau.grid = c(0.1, 0.5, 0.9)

L = 5
#M = 10
D = ncol(X)
phi = phi_generator(L, tau.grid)
M = length(tau.grid)
#tau.grid = tau.grid_generator(L, M)
#phi = phi_generator2(L,M)


maxit = 1000

##### Algorithm 2

a_0 = matrix(1:(D*L), nrow = D)
a_vec_0 = as.vector(a_0)
eta_base = matrix(1, N,M)
eta_0 = eta_vec_gen(y, X, a_vec_0, phi)
a_vec_k = a_vec_0
eta_k = eta_0
R_old = R_vec(a_vec_0, y, X, tau.grid, phi, eta_k)

d_index = vector()
l_index = vector()
for (d in 1:D){
  for (l in 1:L){
    j_equivalent = (l-1)*D + d
    d_index[j_equivalent] = d
    l_index[j_equivalent] = l
    
  }
}

t_start2 = Sys.time()
for(i in 1:maxit){
  #result  = optim(a_vec_k, R_vec2, y=y, X=X, tau.grid=tau.grid, phi=phi, eta=eta_k, d_index=d_index, l_index = l_index, method = "BFGS", gr = gr)
  #a_vec_k = result$par
  a_vec_k = a_vec_opt(y, X, phi, tau.grid, eta_k)
  eta_k = eta_vec_gen(y, X, a_vec_k, phi)
  R_k = R_vec(a_vec_k, y, X, tau.grid, phi, eta_k)
  print(i)
  print(R_old)
  print(R_k)
  if (convergence(R_old, R_k)) break
  R_old = R_k
}
t_end2 = Sys.time()
time_Rvec2 = t_end2 - t_start2

a_res = matrix(a_vec_k, nrow=D)
B = format(a_res%*%phi, scientific = F)

# Compare with CVXR
gLasso = global_qr(taus = tau.grid, phi = phi, X = X, y = y, lambda = 0, lags = 0, f = "gLasso", w = F)
a_CVXR = gLasso$ahat
B_CVXR = a_CVXR%*%phi

# Compare with QR
X_0 = data$X
qrfit = rq(y ~ X_0, tau = tau.grid)
coef(qrfit)
qrfit_coefs = cbind(t(coef(qrfit)), rep(2,M))

