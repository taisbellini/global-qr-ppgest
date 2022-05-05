#### Local MSE ####
# calculates the MSE of the beta_hat by coefficient and by quantile level 
# with the chosen lambda of each iteration according to criteria (AIC/BIC/EBIC)

reps_vector = 1:nrep
clean_reps = reps_vector[! reps_vector %in% error_reps]
nrep_clean = length(clean_reps)
if(nrep_clean >= 200){
  clean_reps = clean_reps[1:200]
  nrep_clean = 200
} 
lambdas_clean = 30

mse_D = array(0, dim = c(D,length(methods_clean), length(criterias)))
mse_M = array(0, dim = c(M,length(methods_clean), length(criterias)))

for (method in methods){
  index = which(methods_clean== method)
  se_aic = numeric(D)
  se_bic = numeric(D)
  se_ebic = numeric(D)
  
  se_aic_M = numeric(M)
  se_bic_M = numeric(M)
  se_ebic_M = numeric(M)
 
  for(i in clean_reps){
    error_mat_aic = (opt_datasets[[method]][,,i, 1] - beta_real)^2
    se_aic = se_aic + apply(error_mat_aic, 1, sum)
    se_aic_M = se_aic_M + apply(error_mat_aic, 2, sum)
    
    error_mat_bic = (opt_datasets[[method]][,,i, 2] - beta_real)^2
    se_bic = se_bic + apply(error_mat_bic, 1, sum)
    se_bic_M = se_bic_M + apply(error_mat_bic, 2, sum)
    
    error_mat_ebic = (opt_datasets[[method]][,,i, 3] - beta_real)^2
    se_ebic = se_ebic + apply(error_mat_ebic, 1, sum)
    se_ebic_M = se_ebic_M + apply(error_mat_ebic, 2, sum)
  }
  mse_D[,index, 1] = se_aic/nrep_clean
  mse_D[,index, 2] = se_bic/nrep_clean
  mse_D[,index, 3] = se_ebic/nrep_clean
  
  mse_M[,index, 1] = se_aic_M/nrep_clean
  mse_M[,index, 2] = se_bic_M/nrep_clean
  mse_M[,index, 3] = se_ebic_M/nrep_clean
}

## Calculate MSE sum and mean per method
#evaluate non-zero coefficients to avoid distortion from computational errors

#mse_sum: criteria x method with MSE sum across all reps using opt lambda for each rep
mse_sum = array(0, dim = c(length(criterias),length(methods_clean)))
colnames(mse_sum) = methods_clean
rownames(mse_sum) = criterias

for (method in methods_clean){
  index = which(methods_clean == method)
  mse_sum[,index] = apply(mse_D[,index,], 2, sum)
}

print(mse_sum)

#mse_mean: criteria x method with MSE mean across all reps using opt lambda for each rep
mse_mean = array(0, dim = c(length(criterias),length(methods_clean)))
colnames(mse_mean) = methods_clean
rownames(mse_mean) = criterias

for (method in methods_clean){
  index = which(methods_clean == method)
  mse_mean[,index] = apply(mse_D[,index,], 2, mean)
}

print(mse_mean)

#### MSE "global" ####

mse_global = array(0, dim = c(lambdas_clean,D,length(methods_clean)))
mse_global_M = array(0, dim = c(lambdas_clean,M,length(methods_clean)))

for (method in methods_clean){
  index = which(methods_clean == method)
  for (l in 1:lambdas_clean){
    se = numeric(D)
    se_M = numeric(M)
    for(i in clean_reps){
      error_mat = (datasets[[method]][,,l,i] - beta_real)^2
      se = se + apply(error_mat, 1, sum)
      se_M = se_M + apply(error_mat, 2, sum)
      
    }
    mse_global[l, ,index] = se/nrep_clean
    mse_global_M[l, , index] = se_M/nrep_clean
  }
}

## Find lowest MSE sum and mean per method
#evaluate non-zero coefficients to avoid distortion from computational errors

#mse_global_sum: lambda x method with the MSE sum across all reps for each lambda.
mse_global_sum = array(0, dim = c(lambdas_clean,length(methods_clean)))
colnames(mse_global_sum) = methods_clean

# sum the MSE of each lambda
for (method in methods_clean){
  index = which(methods_clean == method)
  mse_global_sum[,index] = apply(mse_global[,,index], 1, sum)
}

#mse_global_mean: lambda x method with the MSE mean across all reps for each lambda.
mse_global_mean = array(0, dim = c(lambdas_clean,length(methods_clean)))
colnames(mse_global_mean) = methods_clean

# mean the MSE of each lambda
for (method in methods_clean){
  index = which(methods_clean == method)
  mse_global_mean[,index] = apply(mse_global[,,index], 1, mean)
}

# Global lambda choice

# using sum to find the lambda that generates the lowest MSE across all reps.
mse_global_min_lambda = apply(mse_global_sum, 2, function(mse_vals){
  lambda_index = which(mse_vals == min(mse_vals))
  return(lambdas_to_try[lambda_index])
})
mse_min_lambda_index = apply(mse_global_sum, 2, which.min)

## build data frame with MSE from BIC lambdas and globally chosen lambda for each method

mse_min_global = array(0, dim = c(length(methods_clean), D))

for (method in methods_clean){
  met = which(methods_clean == method)
  mse_min_global[met,] = mse_global[mse_min_lambda_index[met],,met]
}

# Global MSE sum
mse_min_global_sum = apply(mse_min_global, 1, sum)
names(mse_min_global_sum) = methods_clean

mse_min_sum_mat = cbind(round(mse_sum[2,], 4), round(mse_sum[1,], 4), round(mse_min_global_sum, 4))
colnames(mse_min_sum_mat) = c("BIC MSE", "AIC MSE", "Omni MSE")

# Global MSE mean
mse_min_global_mean = apply(mse_min_global, 1, mean)
names(mse_min_global_mean) = methods_clean

mse_min_mean_mat = cbind(round(mse_mean[2,], 4), round(mse_mean[1,], 4), round(mse_min_global_mean, 4))
colnames(mse_min_mean_mat) = c("BIC MSE", "AIC MSE", "Omni MSE")

fig_name = sprintf("fig_outputs/mse_table_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width=350,height=200,bg = "white")
grid.table(mse_min_mean_mat)
dev.off()

fig_name = sprintf("fig_outputs/mse_sum_table_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width=350,height=200,bg = "white")
grid.table(mse_min_sum_mat)
dev.off()

#### Elicitability ####

#### Local Elicitability (/BIC/)

elic = array(0, dim = c(length(methods_clean), length(criterias)))

for (method in methods_clean){
  index = which(methods_clean == method)
    elic_sum_bic = 0
    elic_sum_aic = 0
    for (m in 1:M){
      rho_sum_bic = 0
      rho_sum_aic = 0
      for(i in clean_reps){
        u_bic = drop(ds_y[1,i] - t(ds_X[1,,i])%*%opt_datasets[[method]][,m,i,2])
        rho_sum_bic  = rho_sum_bic + rho_function(u_bic, tau.grid[m])
        
        u_aic = drop(ds_y[1,i] - t(ds_X[1,,i])%*%opt_datasets[[method]][,m,i,1])
        rho_sum_aic  = rho_sum_aic + rho_function(u_aic, tau.grid[m])
      }
      elic_sum_bic = elic_sum_bic + rho_sum_bic/nrep_clean
      elic_sum_aic = elic_sum_aic + rho_sum_aic/nrep_clean
    }
    elic[index, 2] = elic_sum_bic
    elic[index, 1] = elic_sum_aic
}
rownames(elic) = methods_clean


### Global Elic (choose best of all lambdas) ###

elic_global = array(0, dim = c(length(methods_clean),lambdas_clean))

for (method in methods_clean){
  index = which(methods_clean == method)
  for (l in 1:lambdas_clean){
    elic_sum = 0
    for (m in 1:M){
      rho_sum = 0
      for(i in clean_reps){
        u = drop(ds_y[1,i] - t(ds_X[1,,i])%*%datasets[[method]][,m,l,i])
        rho_sum  = rho_sum + rho_function(u, tau.grid[m])
      }
      elic_sum = elic_sum + rho_sum/nrep_clean
    }
    elic_global[index, l] = elic_sum
  }
}

elic_min = apply(elic_global, 1, min)
names(elic_min) = methods_clean
elic_min_lambda_index = apply(elic_global, 1, which.min)

elic_min_lambda = apply(abs(elic), 1, function(elic_vals){
  lambda_index = which.min(elic_vals)
  return(lambdas_to_try[lambda_index])
})

## build data frame with BIC elic and only the best lambda Elic of each method

elic_best_mat = cbind(round(elic[,2], 4), round(elic[,1], 4), round(elic_min,4))
colnames(elic_best_mat) = c("Elicibility BIC", "Elicibility AIC", "Elicitibility Omni")

print(elic_best_mat)
fig_name = sprintf("fig_outputs/elic_table_N%0.0f_M%0.0f_delta%0.2f_s%0.1f.png", N, M, delta, s)
png(fig_name, width=350,height=200,bg = "white")
grid.table(elic_best_mat)
dev.off()

### Plot the histogram of chosen lambdas from AIC and BIC crteria vs. globally chosen lambda

bic_lambdas_ds = data.frame(t(opt_criterias_lambdas[methods_clean_index,2,clean_reps]))
colnames(bic_lambdas_ds) = methods_clean

aic_lambdas_ds = data.frame(t(opt_criterias_lambdas[methods_clean_index,1,clean_reps]))
colnames(aic_lambdas_ds) = methods_clean

filename = sprintf("fig_outputs/lambda_choice_lasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=lasso)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[1], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[1], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[1], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[1], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/lambda_choice_adaLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=`ada lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[2], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[2], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[2], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[2], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Ada LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/lambda_choice_groupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=`group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[3], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[3], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[3], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[3], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Group LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/lambda_choice_adaGroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=`ada group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[4], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[4], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") + 
  #geom_text(aes(x=mse_global_min_lambda[4], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[4], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Ada group LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/lambda_choice_naivegroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=`naive group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[5], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[5], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[7], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[7], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Naive group LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/lambda_choice_naiveAdaGroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=`naive ada group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[6], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[6], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[8], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[8], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Naive ada group LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

#### AIC #####
filename = sprintf("fig_outputs/lambda_choice_AIC_lasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=lasso)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[1], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[1], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[1], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[1], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/lambda_choice_AIC_adaLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=`ada lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[2], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[2], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[2], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[2], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Ada LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/lambda_choice_AIC_groupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=`group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[3], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[3], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[3], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[3], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Group LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/lambda_choice_AIC_adaGroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=`ada group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[4], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[4], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") + 
  #geom_text(aes(x=mse_global_min_lambda[4], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[4], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Ada group LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/lambda_choice_AIC_naivegroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=`naive group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[5], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[5], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[7], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[7], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Naive group LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/lambda_choice_AIC_naiveAdaGroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=`naive ada group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = mse_global_min_lambda[6], color="MSE"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = elic_min_lambda[6], color="Elic."), linetype="dashed", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[8], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[8], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Naive ada group LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

#### Plot estimation of coefficients BIC ####
for(d in c(1, 5, 15, 25)){
  for(method in methods_clean){
    mean_values = apply(opt_datasets[[method]][d,,clean_reps[1:50],2],1, mean)
    png(sprintf("fig_outputs/50reps_BIC_%s_N%0.0f_M%0.0f_D%0.0f.png", method, N, M, d), width=1100, height=570)
    y_inf = min(opt_datasets[[method]][d,,,2])
    y_sup = max(opt_datasets[[method]][d,,,2])
    plot(tau.grid, log(1+beta_real[d,]), type = 'l', col='black', lwd=3, lty='dotted', 
         ylab = "Estimated beta", xlab = "Quantile level", main = sprintf("%s N=%0.0f M=%0.0f D=%0.0f BIC", method, N, M, d), ylim = c(y_inf, y_sup),
         cex.main=1.5)
    for(i in clean_reps[1:50]){
      lines(tau.grid, opt_datasets[[method]][d,,i,2], lwd=2, col = rgb(1,0,0,.1))
    }
    lines(tau.grid, mean_values, type = 'l', lty='dotted', lwd=2, col = rgb(0,0,1))
    dev.off()
  }
}

#### Plot estimation of coefficients AIC ####
for(d in c(1, 5, 15, 25)){
  for(method in methods_clean){
    mean_values = apply(opt_datasets[[method]][d,,clean_reps[1:50],1],1, mean)
    png(sprintf("fig_outputs/50reps_AIC_%s_N%0.0f_M%0.0f_D%0.0f.png", method, N, M, d), width=1100, height=570)
    y_inf = min(opt_datasets[[method]][d,,,1])
    y_sup = max(opt_datasets[[method]][d,,,1])
    plot(tau.grid, log(1+beta_real[d,]), type = 'l', col='black', lwd=3, lty='dotted', 
         ylab = "Estimated beta", xlab = "Quantile level", main = sprintf("%s N=%0.0f M=%0.0f D=%0.0f AIC", method, N, M, d), ylim = c(y_inf, y_sup),
         cex.main=1.5)
    for(i in clean_reps[1:50]){
      lines(tau.grid, opt_datasets[[method]][d,,i,1], lwd=2, col = rgb(1,0,0,.1))
    }
    lines(tau.grid, mean_values, type = 'l', lty='dotted', lwd=2, col = rgb(0,0,1))
    dev.off()
  }
}

### Plot estimation of coefficients Omni ##
for(d in c(1, 5, 15, 25)){
  for(method in methods_clean){
    lambda_opt = mse_min_lambda_index[method]
    mean_values = apply(datasets[[method]][d,,lambda_opt,clean_reps[1:50]],1, mean)
    png(sprintf("fig_outputs/50reps_Omni_%s_N%0.0f_M%0.0f_D%0.0f.png", method, N, M, d))
    y_inf = min(datasets[[method]][d,,lambda_opt,])
    y_sup = max(datasets[[method]][d,,lambda_opt,])
    plot(tau.grid, beta_real[d,], type = 'l', col='black', lwd=3, lty='dotted', 
         ylab = "Estimated beta", xlab = "Quantile level", main = sprintf("%s N=%0.0f M=%0.0f D=%0.0f Omni MSE", method, N, M, d), ylim = c(y_inf, y_sup), 
         cex.main=1.5)
    for(i in clean_reps[1:50]){
      lines(tau.grid, datasets[[method]][d,,lambda_opt,i], lwd=2, col = rgb(1,0,0,.05))
    }
    lines(tau.grid, mean_values, type = 'l', lty='dotted', lwd=2, col = rgb(0,0,1))
    dev.off()
  }
  
}
