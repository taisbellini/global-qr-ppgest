######################################################################################
# This code generates the following datasets: 
# datasets: raw beta hat from simulations for all lambdas and all reps
#   [[method]] D x M x lambda x nrep
# datasets_ahat: raw ahat from simulations for all lambdas and all reps
#   [[method]] D x L x lambda x nrep
# opt_dataset: beta_hat from opt lambda of each iteration according to AIC/BIC/EBIC
#   [[method]] D x M x nrep x criteria
# opt_criterias_lambdas: chosen lambda of each iteration according to AIC/BIC/EBIC
#   method x criteria x nrep
######################################################################################

###### Variables def ####
D_signal=20
D_strong=9
tol = 1e-4

tau.grid = tau.grid_generator(L, M)
phi = phi_generator2(L,tau.grid)
lambdas_to_try = round(10^seq(-3, 3, length.out = lambdas), 4)
beta_real = t(beta(tau.grid, D, D_signal, D_strong))

### Beta real plot ### 
# Define M = 20 to generate this data ##
#for(d in 1:D_signal){
#  filename = sprintf("../dissertation/templateDissertation/fig/beta_real_D%0.0f.png", d)
#  png(filename)
#  plot(tau.grid, beta_real[d,], main = sprintf("Beta coefficients D = %0.0f", d), pch = 16, xlab = "Quantile level", ylab = sprintf("Beta coefficient D = %0.0f", d))
#  dev.off()
#}

### Create data structure to save results for each method ####

## Get beta estimation and zero based on tol

ds_lasso_betas = ds_lasso
ds_lassoW_betas = ds_lassoW
ds_gLasso_betas = ds_gLasso
ds_gLassoW_betas = ds_gLassoW

for(i in 1:nrep){
  for(l in 1:lambdas){
    ds_lasso[,,l,i][abs(ds_lasso[,,l,i])<tol] = 0
    ds_lassoW[,,l,i][abs(ds_lassoW[,,l,i])<tol] = 0
    ds_gLasso[,,l,i][abs(ds_gLasso[,,l,i])<tol] = 0
    ds_gLassoW[,,l,i][abs(ds_gLassoW[,,l,i])<tol] = 0
    ds_naivelasso[,,l,i][abs(ds_naivelasso[,,l,i])<tol] = 0
    ds_naivelassoW[,,l,i][abs(ds_naivelassoW[,,l,i])<tol] = 0
    ds_naivegLasso[,,l,i][abs(ds_naivegLasso[,,l,i])<tol] = 0
    ds_naivegLassoW[,,l,i][abs(ds_naivegLassoW[,,l,i])<tol] = 0
    
    ds_lasso_betas[,,l,i] = ds_lasso[,,l,i]%*%phi
    ds_lassoW_betas[,,l,i] = ds_lassoW[,,l,i]%*%phi
    ds_gLasso_betas[,,l,i] =  ds_gLasso[,,l,i]%*%phi
    ds_gLassoW_betas[,,l,i] = ds_gLassoW[,,l,i]%*%phi
  }
}

methods = c("lasso", "ada lasso", "group lasso", "ada group lasso", "naive lasso", "naive ada lasso",  "naive group lasso", "naive ada group lasso")
datasets = list(
  "lasso" = ds_lasso_betas, 
  "ada lasso" = ds_lassoW_betas,
  "group lasso" = ds_gLasso_betas,
  "ada group lasso" = ds_gLassoW_betas,
  "naive lasso" = ds_naivelasso,
  "naive ada lasso" = ds_naivelassoW,
  "naive group lasso" = ds_naivegLasso,
  "naive ada group lasso" = ds_naivegLassoW
  
)

datasets_ahat = list(
  "lasso" = ds_lasso, 
  "ada lasso" = ds_lassoW,
  "group lasso" = ds_gLasso,
  "ada group lasso" = ds_gLassoW,
  "naive lasso" = ds_naivelasso,
  "naive ada lasso" = ds_naivelassoW,
  "naive group lasso" = ds_naivegLasso,
  "naive ada group lasso" = ds_naivegLassoW
  
)

#### Opt lambda from AIC/BIC for each it ####

criterias = c("aic", "bic", "ebic")

criteria_ds_lasso = array(0, dim=c(nrep, lambdas, length(criterias)))
criteria_ds_lassoW = array(0, dim=c(nrep, lambdas, length(criterias)))
criteria_ds_gLasso = array(0, dim=c(nrep, lambdas, length(criterias)))
criteria_ds_gLassoW = array(0, dim=c(nrep, lambdas, length(criterias)))
criteria_ds_naivelasso = array(0, dim=c(nrep, lambdas, length(criterias)))
criteria_ds_naivelassoW = array(0, dim=c(nrep, lambdas, length(criterias)))
criteria_ds_naivegLasso = array(0, dim=c(nrep, lambdas, length(criterias)))
criteria_ds_naivegLassoW = array(0, dim=c(nrep, lambdas, length(criterias)))

criteria_datasets = list(
  "lasso" = criteria_ds_lasso, 
  "ada lasso" = criteria_ds_lassoW,
  "group lasso" = criteria_ds_gLasso,
  "ada group lasso" = criteria_ds_gLassoW,
  "naive lasso" = criteria_ds_naivelasso,
  "naive ada lasso" = criteria_ds_naivelassoW,
  "naive group lasso" = criteria_ds_naivegLasso,
  "naive ada group lasso" = criteria_ds_naivegLassoW
  
)

## Calculates the dev element equivalent to 3.4 in Frumento for a specific iteration and lambda ##
dev_gen = function(lambda_index, ahat, rep_index, f, w, type){
  y = ds_y[1:N,rep_index]
  X = ds_X[1:N,,rep_index]
  Y = matrix(rep(y,M),N,M)
  TAUS = matrix(rep(tau.grid,N),N,M, byrow = TRUE)
  we = weights_function(y, X, tau.grid, phi, f, w, type)
  l_pen = R_function(ahat, Y, X, TAUS, phi, type = type) + lambdas_to_try[lambda_index]*penalty_function(ahat, f, w, s = 2)
  return(l_pen)
}


## AIC/BIC equivalent to 3.6 and 3.7 in Sotille ##

for(method in methods){
  index = which(methods == method)
  for (i in 1:nrep){
    for(l in 1:lambdas){
      if (is.element(method, c("lasso", "ada lasso", "naive lasso", "naive ada lasso"))){
        f = "lasso"
      }else{
        f = "gLasso"
      }
      if (is.element(method, c("ada lasso", "ada group lasso", "naive ada lasso", "naive ada group lasso"))){
        w = T
      }else{
        w = F
      }
      if (is.element(method, c("lasso", "ada lasso", "group lasso", "ada group lasso"))){
        type = "par"
      }else{
        type = "naive"
      }
      
      ahat = datasets_ahat[[method]][,,l,i]
      dev = dev_gen(l, ahat, i, f, w, type)
      
      mod_coefs = datasets[[method]][,,l,i]
      df = nrow(mod_coefs[rowSums(abs(mod_coefs)) > 0,])
      if (is.null(df)) df = 1
      
      aic = log(dev) + (N^(-1))*2*df
      bic = log(dev) + (N^(-1))*log(N)*df
      ebic = log(dev) + (N^(-1))*log(N)*df*L + (1/2)*df*L*(log(D)/N)
      
      criteria_datasets[[method]][i, l, 1] = aic
      criteria_datasets[[method]][i, l, 2] = bic
      criteria_datasets[[method]][i, l, 3] = ebic
    }
  }
}

### Get betas from lowest AIC, BIC, EBIC ##

opt_ds_lasso = array(0, dim=c(D,L,nrep, length(criterias)))
opt_ds_lassoW = array(0, dim=c(D,L,nrep, length(criterias)))
opt_ds_gLasso = array(0, dim=c(D,L,nrep, length(criterias)))
opt_ds_gLassoW = array(0, dim=c(D,L,nrep, length(criterias)))
opt_ds_naivelasso = array(0, dim=c(D,L,nrep, length(criterias)))
opt_ds_naivelassoW = array(0, dim=c(D,L,nrep, length(criterias)))
opt_ds_naivegLasso = array(0, dim=c(D,L,nrep, length(criterias)))
opt_ds_naivegLassoW = array(0, dim=c(D,L,nrep, length(criterias)))

opt_criterias_lambdas = array(0, dim=c(length(methods), length(criterias), nrep))

opt_datasets = list(
  "lasso" = opt_ds_lasso, 
  "ada lasso" = opt_ds_lassoW,
  "group lasso" = opt_ds_gLasso,
  "ada group lasso" = opt_ds_gLassoW,
  "naive lasso" = opt_ds_naivelasso,
  "naive ada lasso" = opt_ds_naivelassoW,
  "naive group lasso" = opt_ds_naivegLasso,
  "naive ada group lasso" = opt_ds_naivegLassoW
  
)

### Get betas from lowest AIC, BIC, EBIC - no clean filter ##

all_opt_ds_lasso = array(0, dim=c(D,L,nrep, length(criterias)))
all_opt_ds_lassoW = array(0, dim=c(D,L,nrep, length(criterias)))
all_opt_ds_gLasso = array(0, dim=c(D,L,nrep, length(criterias)))
all_opt_ds_gLassoW = array(0, dim=c(D,L,nrep, length(criterias)))
all_opt_ds_naivelasso = array(0, dim=c(D,L,nrep, length(criterias)))
all_opt_ds_naivelassoW = array(0, dim=c(D,L,nrep, length(criterias)))
all_opt_ds_naivegLasso = array(0, dim=c(D,L,nrep, length(criterias)))
all_opt_ds_naivegLassoW = array(0, dim=c(D,L,nrep, length(criterias)))

all_opt_criterias_lambdas = array(0, dim=c(length(methods), length(criterias), nrep))

all_opt_datasets = list(
  "lasso" = opt_ds_lasso, 
  "ada lasso" = opt_ds_lassoW,
  "group lasso" = opt_ds_gLasso,
  "ada group lasso" = opt_ds_gLassoW,
  "naive lasso" = opt_ds_naivelasso,
  "naive ada lasso" = opt_ds_naivelassoW,
  "naive group lasso" = opt_ds_naivegLasso,
  "naive ada group lasso" = opt_ds_naivegLassoW
  
)

for (method in methods){
  index = which(methods == method)
  for (i in 1:nrep){
    # get index from higher lambda that gives lowest aic/bic/ebic values
    aic_opt_lambda_index = max(which(criteria_datasets[[method]][i,,1] == min(criteria_datasets[[method]][i,1:lambdas_clean,1])))
    bic_opt_lambda_index = max(which(criteria_datasets[[method]][i,,2] == min(criteria_datasets[[method]][i,1:lambdas_clean,2])))
    ebic_opt_lambda_index = max(which(criteria_datasets[[method]][i,,3] == min(criteria_datasets[[method]][i,1:lambdas_clean,3])))
    
    # save the corresponding betahat value for opt lambda of each iteration
    opt_datasets[[method]][,,i,1] = datasets[[method]][,,aic_opt_lambda_index,i]
    opt_datasets[[method]][,,i,2] = datasets[[method]][,,bic_opt_lambda_index,i]
    opt_datasets[[method]][,,i,3] = datasets[[method]][,,ebic_opt_lambda_index,i]
    
    # save the chosen lambda for each iteraction
    opt_criterias_lambdas[index, 1, i] = lambdas_to_try[aic_opt_lambda_index]
    opt_criterias_lambdas[index, 2, i] = lambdas_to_try[bic_opt_lambda_index]
    opt_criterias_lambdas[index, 3, i] = lambdas_to_try[ebic_opt_lambda_index]
    
    ## No clean filter ##
    # get index from higher lambda that gives lowest aic/bic/ebic values
    all_aic_opt_lambda_index = max(which(criteria_datasets[[method]][i,,1] == min(criteria_datasets[[method]][i,,1])))
    all_bic_opt_lambda_index = max(which(criteria_datasets[[method]][i,,2] == min(criteria_datasets[[method]][i,,2])))
    all_ebic_opt_lambda_index = max(which(criteria_datasets[[method]][i,,3] == min(criteria_datasets[[method]][i,,3])))
    
    # save the corresponding betahat value for opt lambda of each iteration
    all_opt_datasets[[method]][,,i,1] = datasets[[method]][,,all_aic_opt_lambda_index,i]
    all_opt_datasets[[method]][,,i,2] = datasets[[method]][,,all_bic_opt_lambda_index,i]
    all_opt_datasets[[method]][,,i,3] = datasets[[method]][,,all_ebic_opt_lambda_index,i]
    
    # save the chosen lambda for each iteraction
    all_opt_criterias_lambdas[index, 1, i] = lambdas_to_try[all_aic_opt_lambda_index]
    all_opt_criterias_lambdas[index, 2, i] = lambdas_to_try[all_bic_opt_lambda_index]
    all_opt_criterias_lambdas[index, 3, i] = lambdas_to_try[all_ebic_opt_lambda_index]
  }
}

