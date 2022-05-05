
data_filename = "data_N1000_nrep1000_official.RData"
load(file=data_filename)

###### Variables def ####
#nrep = 10000
#lambdas = 50
if(by10){
  nrep = nrep_total + 10
}else {
  nrep = nrep_total
}
lambdas = 50
D = 30

##### Create data structures ####

ds_lasso = array(0, dim=c(D,L, lambdas, nrep))
ds_lassoW = array(0, dim=c(D,L, lambdas, nrep))
ds_gLasso = array(0, dim=c(D,L, lambdas, nrep))
ds_gLassoW = array(0, dim=c(D,L, lambdas, nrep))
ds_naivelasso = array(0, dim=c(D,L, lambdas, nrep))
ds_naivelassoW = array(0, dim=c(D,L, lambdas, nrep))
ds_naivegLasso = array(0, dim=c(D,L, lambdas, nrep))
ds_naivegLassoW = array(0, dim=c(D,L, lambdas, nrep))
ds_y = Y_array[1:N, 1:nrep]
ds_X = X_array[1:N, , 1:nrep]

for(nr in nrep_s_vec){
  
  nrep_s = nr
  nrep_e = (nr+by)-1
  filename = sprintf("results/results_N%0.0f_M%0.0f_delta%0.2f_s%0.1f_from%0.0f_to%0.0f.RData", N, M, delta, s, nrep_s, nrep_e)
  load(file=filename)
  
  for(i in nrep_s:nrep_e){
    ds_lasso[,,,i] = array(sim[,colnames(sim) == "lasso_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
    ds_lassoW[,,,i] = array(sim[,colnames(sim) == "lassoW_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
    ds_gLasso[,,,i] = array(sim[,colnames(sim) == "gLasso_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
    ds_gLassoW[,,,i] = array(sim[,colnames(sim) == "gLassoW_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
    ds_naivelasso[,,,i] = array(sim[,colnames(sim) == "naivelasso_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
    ds_naivegLasso[,,,i] = array(sim[,colnames(sim) == "naivegLasso_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
    ds_naivelassoW[,,,i] = array(sim[,colnames(sim) == "naivelassoW_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
    ds_naivegLassoW[,,,i] = array(sim[,colnames(sim) == "naivegLassoW_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
  }
  
}


### By 10 #####

if (by10){
  for(nr in nrep_s_vec10){
    
    nrep_s = nr
    nrep_e = (nr+10)-1
    filename = sprintf("results/results_N%0.0f_M%0.0f_delta%0.2f_s%0.1f_from%0.0f_to%0.0f.RData", N, M, delta, s, nrep_s, nrep_e)
    load(file=filename)
    
    for(i in nrep_s:nrep_e){
      ds_lasso[,,,i] = array(sim[,colnames(sim) == "lasso_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
      ds_lassoW[,,,i] = array(sim[,colnames(sim) == "lassoW_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
      ds_gLasso[,,,i] = array(sim[,colnames(sim) == "gLasso_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
      ds_gLassoW[,,,i] = array(sim[,colnames(sim) == "gLassoW_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
      ds_naivelasso[,,,i] = array(sim[,colnames(sim) == "naivelasso_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
      ds_naivegLasso[,,,i] = array(sim[,colnames(sim) == "naivegLasso_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
      ds_naivelassoW[,,,i] = array(sim[,colnames(sim) == "naivelassoW_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
      ds_naivegLassoW[,,,i] = array(sim[,colnames(sim) == "naivegLassoW_CVXR"][,(i-nrep_s)+1], c(D, L, lambdas))
    }
    
  }
  
}



