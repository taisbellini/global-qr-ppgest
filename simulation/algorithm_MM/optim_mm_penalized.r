#### data generation

source("../core/data_generator.R")
source("../core/utils.R")
source("../core/global_qr_CVXR.R")
source("global_qr.R")
### Simulate cross section Data ####
N = 1000
Dstar = 15 # vars != 0
D = 20 # total vars

# data gen function from data_generator.R
data = data_gen(N=N, Dstar=Dstar, D=D)

#beta_vec = c(rep(2, 10), rep(0, 20), rep(1/400, 10))
#data = data_gen_var(epsilon, n = n, beta = beta_vec)

X = data$X
X = cbind(1,X)
y = data$y

# params def
#tau.grid = seq(from = .1, to = .9, by = .1)
#tau.grid = c(0.1, 0.5, 0.9)

D = ncol(X)
#L = 5
#phi = phi_generator(L, tau.grid)
L = 10
M = 10
tau.grid = tau.grid_generator(L, M)
phi = phi_generator2(L,tau.grid)

##### LASSO ####

a_0 = matrix(1:(D*L), nrow = D)

t_start2 = Sys.time()
result = global_qr(a_0, tau.grid = tau.grid, phi = phi, X = X, y = y, lambda = 100, f = "lasso", w = F, criteria = "bic")
t_end2 = Sys.time()
time_Rvec2 = t_end2 - t_start2


a_res = result$ahat_tol
B = result$bhat

round(B, 3)
round(a_res, 3)

# Compare with CVXR
t_start_CVXR = Sys.time()
lasso_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = "opt_v2", f = "lasso", w = T, criteria = "bic")
t_end_CVXR = Sys.time()
time_CVXR = t_end_CVXR - t_start_CVXR
a_CVXR = lasso_CVXR$ahat
B_CVXR = a_CVXR%*%phi

round(B_CVXR, 3)

# Compare with QR
X_0 = data$X
qrfit = rq(y ~ X_0, tau = tau.grid)
B_QR = coef(qrfit)

round(a_res - a_CVXR, 5)
round(B - B_CVXR, 5)
round(B - B_QR, 5)


#### TEST LAMBDA BEHAVIOR LASSO####


a_0 = matrix(1:(D*L), nrow = D)
lambdas_to_try = c(.001, .1, .5, 1, 1.5, 2, 5, 10, 20, 30, 50, 100)

B_results = list()
a_results = list()
abs_a = c()

for(lambda in lambdas_to_try){
  result = global_qr(a_0, tau.grid = tau.grid, phi = phi, X = X, y = y, lambda = lambda, f = "lasso", w = F, criteria = "bic")
  
  a_res = result$ahat_tol
  B = result$bhat
  
  for(i in 1:ncol(X)){
    if(lambda == lambdas_to_try[1]){
      a_results[[i]] = a_res[i,]
      B_results[[i]] = B[i,]
    } else {
      a_results[[i]] = rbind(a_results[[i]], a_res[i,])
      B_results[[i]] = rbind(B_results[[i]], B[i,])
    }
  }
  abs_a = c(abs_a, sum(abs(a_res)))
}

# CVXR

B_results_CVXR = list()
a_results_CVXR = list()
abs_a_CVXR = c()

for(lambda in lambdas_to_try){
  lasso_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = lambda, f = "lasso", w = F)
  
  a_CVXR = lasso_CVXR$ahat
  B_CVXR = a_CVXR%*%phi
  
  for(i in 1:ncol(X)){
    if(lambda == lambdas_to_try[1]){
      a_results_CVXR[[i]] = a_CVXR[i,]
      a_results_CVXR[[i]] = B_CVXR[i,]
    } else {
      a_results_CVXR[[i]] = rbind(a_results_CVXR[[i]], a_CVXR[i,])
      a_results_CVXR[[i]] = rbind(a_results_CVXR[[i]], B_CVXR[i,])
    }
  }
  abs_a_CVXR = c(abs_a_CVXR, sum(abs(a_CVXR)))
}

plot(lambdas_to_try, abs_a)
plot(lambdas_to_try, abs_a_CVXR)

##### GROUP LASSO ####

a_0 = matrix(1:(D*L), nrow = D)

t_start2 = Sys.time()
result = global_qr(a_0, tau.grid = tau.grid, phi = phi, X = X, y = y, lambda = 20, f = "gLasso", w = F, criteria = "bic")
t_end2 = Sys.time()
time_Rvec2 = t_end2 - t_start2

a_res = result$ahat_tol
B = result$bhat

# Compare with CVXR
lasso_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = 20, lags = 0, f = "gLasso", w = F)
a_CVXR = lasso_CVXR$ahat
B_CVXR = a_CVXR%*%phi

# Compare with QR
X_0 = data$X
qrfit = rq(y ~ X_0, tau = tau.grid)
B_QR = coef(qrfit)

round(B - B_CVXR, 3)
round(B - B_QR, 3)

