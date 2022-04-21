#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Inform parameters: N, M, delta, s, nrep_s, nrep_e, lambdas", call.=FALSE)
} else if (length(args)>0) {
  print(args)
}

N = as.numeric(args[1])
M = L = as.numeric(args[2])
nrep_s = as.numeric(args[5])
nrep_e = as.numeric(args[6])
lambdas = ifelse(is.na(args[7]), 100, as.numeric(args[7]))
lasso_lambda = "opt"
gLasso_lambda = "opt"

#load libs

packages <- c("quantreg", "CVXR", "doParallel")

install.packages(setdiff(packages, rownames(installed.packages())), repos = "https://cran-r.c3sl.ufpr.br/")

library(quantreg)
library(CVXR)
library(foreach)
library(doParallel)

data_filename = "data_N1000_nrep1000_official.RData"
#data_filename = "dataStrong_N1000_nrep50.RData"
load(file=data_filename)

fun_vec = c('global_qr_CVXR_opt', 'R_CVXR', 'weights_CVXR', 'penalty_CVXR')
pack_vec = c("quantreg", "doParallel", "CVXR")

# Load aux files
source("core/data_generator.R")
source("core/utils.R")
source("core/global_qr_CVXR.R")

## Opt params
tau.grid = tau.grid_generator(L, M, delta = delta)
phi = phi_generator2(L,tau.grid)
lambdas_to_try = round(10^seq(-3, 3, length.out = lambdas), 4)

# Monte Carlo

start = Sys.time()

numCores=round((detectCores()/4)*3)
cl = makeCluster(numCores)
registerDoParallel(cl)

#for (i in 1:nrep){
sim = foreach(i=nrep_s:nrep_e,  .export= fun_vec, .packages = pack_vec, .combine = "cbind") %dopar% {
  
  X = X_array[1:N,,i]
  y = Y_array[1:N,i]

  lasso_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = lasso_lambda, f = "lasso", s=s, w = F, lambdas = lambdas_to_try, type = "par")
  
  gLasso_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = gLasso_lambda, f = "gLasso", s=s, w = F, lambdas = lambdas_to_try, type = "par")
  
  lassoW_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = lasso_lambda, f = "lasso", s=s, w = T, lambdas = lambdas_to_try, type = "par")
  
  gLassoW_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = gLasso_lambda, f = "gLasso", s=s, w = T, lambdas = lambdas_to_try, type = "par")
  
  naivelasso_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = lasso_lambda, f = "lasso", s=s, w = F, lambdas = lambdas_to_try, type = "naive")
  
  naivegLasso_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = gLasso_lambda, f = "gLasso", s=s, w = F, lambdas = lambdas_to_try, type = "naive")
  
  naivelassoW_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = lasso_lambda, f = "lasso", s=s, w = T, lambdas = lambdas_to_try, type = "naive")
  
  naivegLassoW_CVXR = global_qr_CVXR(taus = tau.grid, phi = phi, X = X, y = y, lambda = gLasso_lambda, f = "gLasso", s=s, w = T, lambdas = lambdas_to_try, type = "naive")

  return(cbind(lasso_CVXR, gLasso_CVXR, lassoW_CVXR, gLassoW_CVXR, naivelasso_CVXR, naivegLasso_CVXR, naivelassoW_CVXR, naivegLassoW_CVXR))
  
}

stopCluster(cl)
  
end = Sys.time()
time = end - start
print(time)

filename = sprintf("results/results_N%0.0f_M%0.0f_from%0.0f_to%0.0f.RData", N, M, delta, s, nrep_s, nrep_e)
#filename = sprintf("results/resultsbetaStrong_N%0.0f_M%0.0f_from%0.0f_to%0.0f.RData", N, M, delta, s, nrep_s, nrep_e)
save(sim, file = filename)





