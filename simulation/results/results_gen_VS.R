reps_vector = 1:nrep
clean_reps = reps_vector[! reps_vector %in% error_reps]
nrep_clean = length(clean_reps)
if(nrep_clean >= 200){
  clean_reps = clean_reps[1:200]
  nrep_clean = 200
} 
lambdas_clean = 30

#### FVCI ####
#Average fraction of relevant covariates included and irrelevant covariates excluded;

### Local BIC ##########

FVCI = array(0, dim = c(length(methods_clean), length(criterias)))
rownames(FVCI) = methods_clean

for (method in methods_clean){
  index = which(methods_clean == method)
    fvci_reps_BIC = numeric(nrep)
    fvci_reps_AIC = numeric(nrep)
    
    for(i in clean_reps){
      d_in_right_BIC = 0
      d_out_right_BIC = 0
      
      d_in_right_AIC = 0
      d_out_right_AIC = 0

      for(d_in in 1:D_signal){
        if (sum(abs(opt_datasets[[method]][d_in,,i,2])> 0) > 0) d_in_right_BIC = d_in_right_BIC + 1
        if (sum(abs(opt_datasets[[method]][d_in,,i,1]) > 0) > 0) d_in_right_AIC = d_in_right_AIC + 1
        
      }
      for(d_out in (D_signal+1):D){
        if (sum(abs(opt_datasets[[method]][d_out,,i,2]) < tol) == M) d_out_right_BIC = d_out_right_BIC + 1
        if (sum(abs(opt_datasets[[method]][d_out,,i,1]) < tol) == M) d_out_right_AIC = d_out_right_AIC + 1
        
      }
      fvci_reps_BIC[i] = (d_in_right_BIC + d_out_right_BIC) / D
      fvci_reps_AIC[i] = (d_in_right_AIC + d_out_right_AIC) / D
      
    }
    fvci_reps_BIC = fvci_reps_BIC[clean_reps]
    fvci_reps_AIC = fvci_reps_AIC[clean_reps]
    FVCI[index, 2] = mean(fvci_reps_BIC)
    FVCI[index, 1] = mean(fvci_reps_AIC)
    
}

##### Omni ##################

FVCI_omni = array(0, dim = c(length(methods_clean),lambdas_clean))

for (method in methods_clean){
  index = which(methods_clean == method)
  for (l in 1:lambdas_clean){
    fvci_reps = numeric(nrep)
    for(i in clean_reps){
      d_in_right = 0
      d_out_right = 0
      for(d_in in 1:D_signal){
        if (sum(abs(datasets[[method]][d_in,,l,i]) > 0) > 0) d_in_right = d_in_right + 1
      }
      for(d_out in (D_signal+1):D){
        if (sum(abs(datasets[[method]][d_out,,l,i])< tol) == M) d_out_right = d_out_right + 1
      }
      fvci_reps[i] = (d_in_right + d_out_right) / D
    }
    fvci_reps = fvci_reps[clean_reps]
    FVCI_omni[index, l] = mean(fvci_reps)
  }
}

FVCI_max = apply(FVCI_omni, 1, max)
names(FVCI_max) = methods_clean

FVCI_max_lambda = apply(FVCI_omni, 1, function(FVCI_vals){
  lambda_index = max(which(FVCI_vals == max(FVCI_vals)))
  return(lambdas_to_try[lambda_index])
})
names(FVCI_max_lambda) = methods_clean

FVCI_max_mat = cbind(round(FVCI[,2],4), round(FVCI[,1],4), round(FVCI_max,4))
colnames(FVCI_max_mat) = c("FVCI BIC", "FVCI AIC", "FVCI Omni")

fig_name = sprintf("fig_outputs/FVCI_table_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width=350,height=200,bg = "white")
grid.table(FVCI_max_mat)
dev.off()

####TMI ####
#Fraction of replications where all relevant covariates were included;

#### TMI BIC #####

TMI = array(0, dim = c(length(methods_clean), length(criterias)))
rownames(TMI) = methods_clean

for (method in methods_clean){
  index = which(methods_clean == method)
  TMI_reps_BIC = 0
  TMI_reps_AIC = 0
  
  for(i in clean_reps){
    d_in_right_BIC = 0
    d_in_right_AIC = 0
    
    for(d_in in 1:D_signal){
      if (sum(abs(opt_datasets[[method]][d_in,,i,2])> 0) > 0) d_in_right_BIC = d_in_right_BIC + 1
      if (sum(abs(opt_datasets[[method]][d_in,,i,1])> 0) > 0) d_in_right_AIC = d_in_right_AIC + 1
      
    }
    if(d_in_right_BIC == D_signal) TMI_reps_BIC = TMI_reps_BIC + 1
    if(d_in_right_AIC == D_signal) TMI_reps_AIC = TMI_reps_AIC + 1
    
    TMI[index, 2] = TMI_reps_BIC/nrep_clean
    TMI[index, 1] = TMI_reps_AIC/nrep_clean
    
  }
}


#### TMI Omni #####

TMI_omni = array(0, dim = c(length(methods_clean),lambdas_clean))

for (method in methods_clean){
  index = which(methods_clean == method)
  for (l in 1:lambdas_clean){
    TMI_reps = 0
    for(i in clean_reps){
      d_in_right = 0
      for(d_in in 1:D_signal){
        if (sum(abs(datasets[[method]][d_in,,l,i])> 0) > 0) d_in_right = d_in_right + 1
      }
      if(d_in_right == D_signal) TMI_reps = TMI_reps + 1
    }
    
    TMI_omni[index, l] = TMI_reps/nrep_clean
  }
}

TMI_max = apply(TMI_omni, 1, max)
names(TMI_max) = methods_clean

TMI_max_lambda = apply(TMI_omni, 1, function(TMI_vals){
  lambda_index = which.max(TMI_vals)
  return(lambdas_to_try[lambda_index])
})

TMI_max_mat = cbind(round(TMI[,2],4),round(TMI[,1],4), round(TMI_max,4))
colnames(TMI_max_mat) = c("TMI BIC","TMI AIC", "TMI Omni")

fig_name = sprintf("fig_outputs/TMI_table_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width=350,height=200,bg = "white")
grid.table(TMI_max_mat)
dev.off()

#### FRVI ####
#Average fraction of relevant covariates included;


#### FRVI BIC ######

FRVI = array(0, dim = c(length(methods_clean), length(criterias)))

for (method in methods_clean){
  index = which(methods_clean == method)
    frvi_reps_BIC = numeric(nrep)
    frvi_reps_AIC = numeric(nrep)
    
    for(i in clean_reps){
      d_in_right_BIC = 0
      d_in_right_AIC = 0
      
      for(d_in in 1:D_signal){
        if (sum(abs(opt_datasets[[method]][d_in,,i,2])> 0) > 0) d_in_right_BIC = d_in_right_BIC + 1
        if (sum(abs(opt_datasets[[method]][d_in,,i,1])> 0) > 0) d_in_right_AIC = d_in_right_AIC + 1
      }
      frvi_reps_BIC[i] = d_in_right_BIC / D_signal
      frvi_reps_AIC[i] = d_in_right_AIC / D_signal
      
    }
    frvi_reps_BIC = frvi_reps_BIC[clean_reps]
    frvi_reps_AIC = frvi_reps_AIC[clean_reps]
    
    FRVI[index,2] = mean(frvi_reps_BIC)
    FRVI[index,1] = mean(frvi_reps_AIC)
}

#### FRVI Omni ######

FRVI_omni = array(0, dim = c(length(methods_clean),lambdas_clean))

for (method in methods_clean){
  index = which(methods_clean == method)
  for (l in 1:lambdas_clean){
    frvi_reps = numeric(nrep)
    for(i in clean_reps){
      d_in_right = 0
      for(d_in in 1:D_signal){
        if (sum(abs(datasets[[method]][d_in,,l,i])> 0) > 0) d_in_right = d_in_right + 1
      }
      frvi_reps[i] = d_in_right / D_signal
    }
    frvi_reps = frvi_reps[clean_reps]
    FRVI_omni[index, l] = mean(frvi_reps)
  }
}

FRVI_max = apply(FRVI_omni, 1, max)
names(FRVI_max) = methods_clean

FRVI_max_lambda = apply(FRVI_omni, 1, function(FRVI_vals){
  lambda_index = which.max(FRVI_vals)
  return(lambdas_to_try[lambda_index])
})

FRVI_max_mat = cbind(round(FRVI[,2],4), round(FRVI[,1],4), round(FRVI_max,4))
colnames(FRVI_max_mat) = c("FRVI BIC","FRVI AIC", "FRVI Omni")

fig_name = sprintf("fig_outputs/FRVI_table_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width=350,height=200,bg = "white")
grid.table(FRVI_max_mat)
dev.off()


#### FIVE ####
#Average fraction of irrelevant covariates excluded;

#### FIVE BIC ####


FIVE = array(0, dim = c(length(methods_clean), length(criterias)))

for (method in methods_clean){
  index = which(methods_clean == method)
  five_reps_BIC = numeric(nrep)
  five_reps_AIC = numeric(nrep)
  
  for(i in clean_reps){
    d_out_right_BIC = 0
    d_out_right_AIC = 0
    
    for(d_out in (D_signal + 1):D){
      if (sum(abs(opt_datasets[[method]][d_out,,i,2])< tol) == M) d_out_right_BIC = d_out_right_BIC + 1
      if (sum(abs(opt_datasets[[method]][d_out,,i,1])< tol) == M) d_out_right_AIC = d_out_right_AIC + 1
      
    }
    five_reps_BIC[i] = d_out_right_BIC / (D-D_signal)
    five_reps_AIC[i] = d_out_right_AIC / (D-D_signal)
    
  }
  five_reps_BIC = five_reps_BIC[clean_reps]
  five_reps_AIC = five_reps_AIC[clean_reps]
  
  FIVE[index,2] = mean(five_reps_BIC)
  FIVE[index,1] = mean(five_reps_AIC)
  
}

##### FIVE Omni #####

FIVE_omni = array(0, dim = c(length(methods_clean),lambdas_clean))

for (method in methods_clean){
  index = which(methods_clean == method)
  for (l in 1:lambdas_clean){
    five_reps = numeric(nrep)
    for(i in clean_reps){
      d_out_right = 0
      for(d_out in (D_signal + 1):D){
        if (sum(abs(datasets[[method]][d_out,,l,i])< tol) == M) d_out_right = d_out_right + 1
      }
      five_reps[i] = d_out_right / (D-D_signal)
    }
    five_reps = five_reps[clean_reps]
    FIVE_omni[index, l] = mean(five_reps)
  }
}

FIVE_max = apply(FIVE_omni, 1, max)
names(FIVE_max) = methods_clean

FIVE_max_lambda = apply(FIVE_omni, 1, function(FIVE_vals){
  lambda_index = which.max(FIVE_vals)
  return(lambdas_to_try[lambda_index])
})

FIVE_max_mat = cbind(round(FIVE[,2],4), round(FIVE[,1],4), round(FIVE_max,4))
colnames(FIVE_max_mat) = c("FIVE BIC","FIVE AIC", "FIVE Omni")

fig_name = sprintf("fig_outputs/FIVE_table_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width=350,height=200,bg = "white")
grid.table(FIVE_max_mat)
dev.off()


#### NIV ####
#Average number of included covariates. 

### BIC NIV ###

NIV = array(0, dim = c(length(methods_clean), length(criterias)))
included = array(0, dim=c(D, length(methods_clean), length(criterias)))

for (method in methods_clean){
  index = which(methods_clean == method)
  niv_reps_BIC = numeric(nrep)
  niv_reps_AIC = numeric(nrep)
  for(i in clean_reps){
    d_in_BIC = 0
    d_in_AIC = 0
    for(d in 1:D){
      if (sum(abs(opt_datasets[[method]][d,,i,2])> 0 ) > 0){
        d_in_BIC = d_in_BIC + 1
        included[d, index, 2] = included[d, index, 2] + 1
      } 
      if (sum(abs(opt_datasets[[method]][d,,i,1])> 0 ) > 0){
        d_in_AIC = d_in_AIC + 1
        included[d, index, 1] = included[d, index, 1] + 1
      } 
    }
    
    niv_reps_BIC[i] = d_in_BIC
    niv_reps_AIC[i] = d_in_AIC
  }
  niv_reps_BIC = niv_reps_BIC[clean_reps]
  niv_reps_AIC = niv_reps_BIC[clean_reps]
  NIV[index, 2] = mean(niv_reps_BIC)
  NIV[index, 1] = mean(niv_reps_AIC)
}

NIV_BIC_mat = cbind(NIV, opt_criterias_lambdas[methods_clean_index])
colnames(NIV_BIC_lambdas_mat) = c("NIV", "BIC lambda")
rownames(NIV_BIC_lambdas_mat) = methods_clean

fig_name = sprintf("fig_outputs/FIVE_table_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width=350,height=200,bg = "white")
grid.table(FIVE_max_mat)
dev.off()

### Plot which variables included ###

bic_included = data.frame(included[,,2])
colnames(bic_included) = methods_clean
bic_included$var = 1:30

aic_included = data.frame(included[,,1])
colnames(aic_included) = methods_clean
aic_included$var = 1:30

fig_name = sprintf("fig_outputs/which_vars_included_lasso_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=bic_included, aes(x = var, y=lasso)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("LASSO N=%0.0f M=%0.0f BIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

fig_name = sprintf("fig_outputs/which_vars_included_ada lasso_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=bic_included, aes(x = var, y=bic_included$`ada lasso`)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("Ada LASSO N=%0.0f M=%0.0f BIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

fig_name = sprintf("fig_outputs/which_vars_included_group lasso_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=bic_included, aes(x = var, y=`group lasso`)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("Group LASSO N=%0.0f M=%0.0f BIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

fig_name = sprintf("fig_outputs/which_vars_included_ada group lasso_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=bic_included, aes(x = var, y=`ada group lasso`)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("Ada Group LASSO N=%0.0f M=%0.0f BIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

fig_name = sprintf("fig_outputs/which_vars_included_naive group lasso_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=bic_included, aes(x = var, y=`naive group lasso`)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("Naive Group LASSO N=%0.0f M=%0.0f BIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

fig_name = sprintf("fig_outputs/which_vars_included_naive ada group lasso_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=bic_included, aes(x = var, y=`naive ada group lasso`)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("Naive Ada Group LASSO N=%0.0f M=%0.0f BIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

### AIC ###

fig_name = sprintf("fig_outputs/which_vars_included_lasso_AIC_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=aic_included, aes(x = var, y=lasso)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("LASSO N=%0.0f M=%0.0f AIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

fig_name = sprintf("fig_outputs/which_vars_included_ada lasso_AIC_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=aic_included, aes(x = var, y=bic_included$`ada lasso`)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("Ada LASSO N=%0.0f M=%0.0f AIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

fig_name = sprintf("fig_outputs/which_vars_included_group lasso_AIC_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=aic_included, aes(x = var, y=`group lasso`)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("Group LASSO N=%0.0f M=%0.0f AIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

fig_name = sprintf("fig_outputs/which_vars_included_ada group lasso_AIC_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=aic_included, aes(x = var, y=`ada group lasso`)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("Ada Group LASSO N=%0.0f M=%0.0f AIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

fig_name = sprintf("fig_outputs/which_vars_included_naive group lasso_AIC_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=aic_included, aes(x = var, y=`naive group lasso`)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("Naive Group LASSO N=%0.0f M=%0.0f AIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

fig_name = sprintf("fig_outputs/which_vars_included_naive ada group lasso_AIC_N%0.0f_M%0.0f.png", N, M)
png(fig_name)
ggplot(data=aic_included, aes(x = var, y=`naive ada group lasso`)) +
  geom_bar(stat="identity") + xlab("Variable") + ylab("# reps included") + ggtitle(sprintf("Naive Ada Group LASSO N=%0.0f M=%0.0f AIC", N, M)) +
  theme(text = element_text(family="serif", size = 18))
dev.off()

### NIV Omni ####

NIV_omni = array(0, dim = c(length(methods_clean),lambdas_clean))

for (method in methods_clean){
  index = which(methods_clean == method)
  for (l in 1:lambdas_clean){
    niv_reps = numeric(nrep)
    for(i in clean_reps){
      d_in = 0
      for(d in 1:D){
        if (sum(abs(datasets[[method]][d,,l,i])) > 0) d_in = d_in + 1
      }
      
      niv_reps[i] = d_in
    }
    niv_reps = niv_reps[clean_reps]
    NIV_omni[index, l] = mean(niv_reps)
  }
}
colnames(NIV_omni) = lambdas_to_try[1:30]
rownames(NIV_omni) = methods_clean

### Plot variable included per lambda ###
## dataframe: method x var x lambda x nreps_included

NIV_df = data.frame(matrix(ncol = 4, nrow = 0))
colnames(NIV_df) = c("method", "var", "lambda", "nreps_included")


for(method in methods_clean){
  for(d in 1:D){
    for (l in 1:lambdas_clean){
      count_included = 0
      for (i in clean_reps) {
        if (sum(abs(datasets[[method]][d,,l,i])> 0 ) > 0) count_included = count_included + 1
        
      }
      NIV_df[nrow(NIV_df)+1,] = c(method, d, l, count_included)
    }
  }
}

NIV_df$method = as.factor(NIV_df$method)
NIV_df$var = factor(NIV_df$var, levels = as.character(1:30))
NIV_df$lambda = factor(NIV_df$lambda, levels = as.character(1:lambdas_clean))
NIV_df$nreps_included = as.numeric(NIV_df$nreps_included)

for(method in methods_clean){
    filename = sprintf("fig_outputs/VS_vars_%s_N%0.0f_M%0.0f.png", method, N, M)
    png(filename)
      plot = ggplot(NIV_df[NIV_df$method == method,], aes(x = lambda, y = var)) +
      geom_point(size = 3, aes(colour=nreps_included)) + 
      scale_color_gradient(low = "blue", high = "red", limits = c(0,200)) + 
      geom_hline(yintercept = 20.5) + 
      geom_hline(yintercept = 9.5) +
      labs(colour="# Rep. Included") + 
      xlab("Lambda") + 
      ylab("Covariate") +
      ggtitle(sprintf("%s N=%0.0f M=%0.0f",method, N, M)) + 
      scale_x_discrete(breaks=c(1,15,30),
                       labels=c(lambdas_to_try[1],lambdas_to_try[15], lambdas_to_try[30])) +
      scale_y_discrete(breaks=c(1,15,30))
      print(plot)
    dev.off()
}

### Correct ####
#Fraction of replications where model was correctly selected;

#### Correct BIC #####

Corr = array(0, dim = c(length(methods_clean)))

for (method in methods_clean){
  index = which(methods_clean == method)
  Corr_reps = 0
  for(i in clean_reps){
    d_in_right = 0
    for(d_in in 1:D_signal){
      if (sum(abs(opt_datasets[[method]][d_in,,i,2])>= tol) > 0) d_in_right = d_in_right + 1
    }
    d_out_right = 0
    for(d_out in (D_signal + 1):D){
      if (sum(abs(opt_datasets[[method]][d_out,,i,2]) < tol) > 0) d_out_right = d_out_right + 1
    }
    if((d_in_right == D_signal) && (d_out_right == (D-D_signal))){
      Corr_reps = Corr_reps + 1
    }
  }
  Corr[index] = Corr_reps/nrep_clean
}

Corr = t(matrix(Corr))
colnames(Corr) = methods_clean
row.names(Corr) = c("Correct")

#### Correct Omni #####

Corr_omni = array(0, dim = c(length(methods_clean),lambdas_clean))

for (method in methods_clean){
  index = which(methods_clean == method)
  for (l in 1:lambdas_clean){
    Corr_reps = 0
    for(i in clean_reps){
      d_in_right = 0
      for(d_in in 1:D_signal){
        if (sum(abs(datasets[[method]][d_in,,l,i])>= tol) > 0) d_in_right = d_in_right + 1
      }
      d_out_right = 0
      for(d_out in (D_signal + 1):D){
        if (sum(abs(datasets[[method]][d_out,,l,i])< tol) > 0) d_out_right = d_out_right + 1
      }
      if((d_in_right == D_signal) && (d_out_right == (D-D_signal))){
        Corr_reps = Corr_reps + 1
      }
    }
    
    Corr_omni[index, l] = Corr_reps/nrep_clean
  }
}

colnames(Corr_omni) = 1:lambdas_clean
rownames(Corr_omni) = methods_clean

fig_name = sprintf("fig_outputs/Corr_BIC_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width = 1000)
grid.table(Corr)
dev.off()

fig_name = sprintf("fig_outputs/Corr_omni_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width = 1500)
grid.table(Corr_omni)
dev.off()


### Overfit ####
#Fraction of replications where model added all correct variables + non-relevant;

#### Overfit BIC #####

Over = array(0, dim = c(length(methods_clean)))

for (method in methods_clean){
  index = which(methods_clean == method)
  Over_reps = 0
  for(i in clean_reps){
    d_in_right = 0
    for(d_in in 1:D_signal){
      if (sum(abs(opt_datasets[[method]][d_in,,i,2]) >= tol) > 0) d_in_right = d_in_right + 1
    }
    d_out_right = 0
    for(d_out in (D_signal + 1):D){
      if (sum(abs(opt_datasets[[method]][d_out,,i,2])< tol) >0) d_out_right = d_out_right + 1
    }
    if((d_in_right == D_signal) && (d_out_right < (D-D_signal))){
      Over_reps = Over_reps + 1
    }
  }
  Over[index] = Over_reps/nrep_clean
}

Over = t(matrix(Over))
colnames(Over) = methods_clean
row.names(Over) = c("Overfit")

#### Overfit Omni #####

Over_omni = array(0, dim = c(length(methods_clean),lambdas_clean))

for (method in methods_clean){
  index = which(methods_clean == method)
  for (l in 1:lambdas_clean){
    Over_reps = 0
    for(i in clean_reps){
      d_in_right = 0
      for(d_in in 1:D_signal){
        if (sum(abs(datasets[[method]][d_in,,l,i])>= tol) > 0) d_in_right = d_in_right + 1
      }
      d_out_right = 0
      for(d_out in (D_signal + 1):D){
        if (sum(abs(datasets[[method]][d_out,,l,i])< tol) > 0) d_out_right = d_out_right + 1
      }
      if((d_in_right == D_signal) && (d_out_right < (D-D_signal))){
        Over_reps = Over_reps + 1
      }
    }
    
    Over_omni[index, l] = Over_reps/nrep_clean
  }
}

colnames(Over_omni) = 1:lambdas_clean
rownames(Over_omni) = methods_clean

fig_name = sprintf("fig_outputs/Over_BIC_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width = 1000)
grid.table(Over)
dev.off()

fig_name = sprintf("fig_outputs/Over_omni_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width = 1500)
grid.table(Over_omni)
dev.off()

### Undefit ####
#Fraction of replications where model added all correct variables + non-relevant;

#### Undefit BIC #####

Under = array(0, dim = c(length(methods_clean)))

for (method in methods_clean){
  index = which(methods_clean == method)
  Under_reps = 0
  for(i in clean_reps){
    d_in_right = 0
    for(d_in in 1:D_signal){
      if (sum(abs(opt_datasets[[method]][d_in,,i,2])>= tol) > 0) d_in_right = d_in_right + 1
    }
    d_out_right = 0
    for(d_out in (D_signal + 1):D){
      if (sum(abs(opt_datasets[[method]][d_out,,i,2])< tol) > 0) d_out_right = d_out_right + 1
    }
    if((d_in_right < D_signal) && (d_out_right == (D-D_signal))){
      Under_reps = Under_reps + 1
    }
  }
  Under[index] = Under_reps/nrep_clean
}

Under = t(matrix(Under))
colnames(Under) = methods_clean
row.names(Under) = c("Underfit")

#### Underfit Omni #####

Under_omni = array(0, dim = c(length(methods_clean),lambdas_clean))

for (method in methods_clean){
  index = which(methods_clean == method)
  for (l in 1:lambdas_clean){
    Under_reps = 0
    for(i in clean_reps){
      d_in_right = 0
      for(d_in in 1:D_signal){
        if (sum(abs(datasets[[method]][d_in,,l,i])>=tol) > 0) d_in_right = d_in_right + 1
      }
      d_out_right = 0
      for(d_out in (D_signal + 1):D){
        if (sum(abs(datasets[[method]][d_out,,l,i])< tol) > 0) d_out_right = d_out_right + 1
      }
      if((d_in_right < D_signal) && (d_out_right == (D-D_signal))){
        Under_reps = Under_reps + 1
      }
    }
    
    Under_omni[index, l] = Under_reps/nrep_clean
  }
}

colnames(Under_omni) = 1:lambdas_clean
rownames(Under_omni) = methods_clean

fig_name = sprintf("fig_outputs/Under_BIC_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width = 1000)
grid.table(Under)
dev.off()

fig_name = sprintf("fig_outputs/Under_omni_N%0.0f_M%0.0f.png", N, M)
png(fig_name, width = 1500)
grid.table(Under_omni)
dev.off()


##### Lambda choice #############

### Plot the histogram of chosen lambdas from BIC crteria vs. globally chosen lambda

bic_lambdas_ds = data.frame(t(opt_criterias_lambdas[,2,]))
aic_lambdas_ds = data.frame(t(opt_criterias_lambdas[,1,]))
colnames(bic_lambdas_ds) = methods
colnames(aic_lambdas_ds) = methods

filename = sprintf("fig_outputs/VS_lambda_choice_lasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=lasso)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[1], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[1], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[1], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[1], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[1], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[1], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/VS_lambda_choice_adaLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=`ada lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[2], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[2], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[2], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[2], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[2], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[2], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Ada LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/VS_lambda_choice_groupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=`group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[3], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[3], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[3], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[3], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[3], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[3], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Group LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/VS_lambda_choice_adaGroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=`ada group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[4], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[4], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[4], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[4], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[4], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[4], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Ada group LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/VS_lambda_choice_naiveGroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=`naive group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[5], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[5], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[5], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[5], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[7], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[7], label="\nElicitability", y=0.2), colour="red", angle=90) +
    ggtitle(sprintf("naive group LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/VS_lambda_choice_naiveAdaGroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
  ggplot(bic_lambdas_ds, aes(x=`naive ada group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[6], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[6], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[6], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[6], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[8], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[8], label="\nElicitability", y=0.2), colour="red", angle=90) +
    ggtitle(sprintf("Naive ada group LASSO N=%0.0f M=%0.0f BIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()


### AIC ###

filename = sprintf("fig_outputs/VS_lambda_choice_AIC_lasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=lasso)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[1], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[1], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[1], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[1], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[1], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[1], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/VS_lambda_choice_AIC_adaLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=`ada lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[2], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[2], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[2], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[2], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[2], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[2], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Ada LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/VS_lambda_choice_AIC_groupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=`group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[3], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[3], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[3], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[3], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[3], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[3], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Group LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/VS_lambda_choice_AIC_adaGroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=`ada group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[4], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[4], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[4], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[4], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[4], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[4], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Ada group LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/VS_lambda_choice_AIC_naiveGroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=`naive group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[5], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[5], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[5], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[5], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[7], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[7], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("naive group LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

filename = sprintf("fig_outputs/VS_lambda_choice_AIC_naiveAdaGroupLasso_N%0.0f_M%0.0f.png", N, M)
png(filename, width = 820, height = 520)
ggplot(aic_lambdas_ds, aes(x=`naive ada group lasso`)) +
  geom_density(fill="snow3") +
  xlim(c(0,4)) +
  #scale_x_break(c(0.3, lambdas_to_try[30]), scales = 1.5) +
  geom_vline(aes(xintercept = FVCI_max_lambda[6], color = "FVCI"), linetype="solid", size = 1) + 
  geom_vline(aes(xintercept = TMI_max_lambda[6], color="TMI"), linetype="dashed", size = 1) + 
  geom_vline(aes(xintercept = FRVI_max_lambda[6], color="FRVI"), linetype="dotted", size = 1) + 
  geom_vline(aes(xintercept = FIVE_max_lambda[6], color="FIVE"), linetype="twodash", size = 1) +
  labs(colour="Omni") +
  #geom_text(aes(x=mse_global_min_lambda[8], label="\nMSE", y=0.2), colour="red", angle=90) +
  #geom_text(aes(x=elic_min_lambda[8], label="\nElicitability", y=0.2), colour="red", angle=90) +
  ggtitle(sprintf("Naive ada group LASSO N=%0.0f M=%0.0f AIC", N, M)) + xlab("lambda") + ylab("Frequency") + 
  theme(text = element_text(family="serif", size = 18))
dev.off()

