### Load datasets already generated ###
source("utils.R")
source("global_qr_CVXR.R")
source("data_generator.R")
library(reshape2)
library(gridExtra)
library(ggplot2)
library(quantreg)
library(CVXR)
library(dplyr)

N = 1000
M = L = 10
delta = 0.05
s = 2
lambdas = 50
#nrep = 600

D = 30
D_signal=20
D_strong=9
tol = 1e-4

tau.grid = tau.grid_generator(L, M)
phi = phi_generator2(L,tau.grid)
lambdas_to_try = round(10^seq(-3, 3, length.out = lambdas), 4)
beta_real = t(beta(tau.grid, D, D_signal, D_strong))

criterias = c("aic", "bic", "ebic")

filename = sprintf("results/datasets_N%0.0f_M%0.0f.RData", N, M)
load(file = filename)

# rename naive to direct

methods = c("lasso", "ada lasso", "group lasso", "ada group lasso", "direct lasso", "direct ada lasso",  "direct group lasso", "direct ada group lasso")
methods_clean = c("lasso", "ada lasso", "group lasso", "ada group lasso", "direct group lasso", "direct ada group lasso")

#names(all_opt_datasets) =
#names(criteria_datasets) =
names(datasets) =
names(datasets_ahat) =
names(opt_datasets) = methods


## generate images with "results_gen_VS.R" and "results_gen_Coefs.R")
