###### Load dependencies

source("../core/data_generator.R")
source("../core/utils.R")
source("../core/global_qr_CVXR.R")
library(reshape2)
library(gridExtra)
library(ggplot2)
library(quantreg)
library(CVXR)
library(dplyr)

### Select the scenario to generate
N = 1000
M = L = 10
delta = 0.05
s = 2

lambdas_clean = 30

### Define the filenames for results_combine
nrep_total = 600
by = 50
nrep_s_vec = seq(from=1, to=nrep_total, by=by)

## Case N500 M10
by10 = F
from_by10 = 351
to_by10 = 360
nrep_s_vec10 = seq(from=from_by10, to=to_by10, by=10)

# Generate datasets: datasets, opt_datasets, criteria_datasets, opt_criterias_lambdas, all_opt_criterias_lambdas  (takes a while)
source("results_combine.R")
source("results_gen.R")

### Plot the lambdas that zero all matrix

# # of reps wwhere that lambda generated an execution error (DxM matrix all zeroes)
error_cases = array(0, dim=c(lambdas, length(methods)))
colnames(error_cases) = methods
rownames(error_cases) = lambdas_to_try

error_cases_BIC = array(0, dim=c(lambdas, length(methods)))
colnames(error_cases_BIC) = methods
rownames(error_cases_BIC) = lambdas_to_try

for (method in methods){
  index = which(methods == method)
  for (l in 1:lambdas){
    for (i in 1:nrep){
      if (sum(abs(datasets[[method]][,,l,i])) == 0) {
        error_cases[l,index] = error_cases[l,index] + 1
        if(l == opt_criterias_lambdas[index,2,i]) {
          error_cases_BIC[l, index] = error_cases_BIC[l, index] + 1
        }
      }
    }
  }
}

error_reps = c()

methods_clean = methods[-c(5,6)]
methods_clean_index = c(1,2,3,4,7,8)

# cleanup replications with errors, see what is left
for(i in 1:nrep){
  stop = FALSE
  for (method in methods_clean){
    index = which(methods == method)
    for (l in 1:lambdas_clean){
      if (sum(abs(datasets[[method]][,,l,i])) == 0) {
        error_reps = c(error_reps, i)
        stop = TRUE
        break
      }
    }
    if(stop) {break}
  }
}
nrep - length(error_reps)

######################################################################################
# Save the datasets to RDataFile: 
# datasets: raw beta hat from simulations for all lambdas and all reps
#   [[method]] D x M x lambda x nrep
# datasets_ahat: raw ahat from simulations for all lambdas and all reps
#   [[method]] D x L x lambda x nrep
# opt_dataset: beta_hat from opt lambda of each iteration according to AIC/BIC/EBIC
#   [[method]] D x M x nrep x criteria
# opt_criterias_lambdas: chosen lambda of each iteration according to AIC/BIC/EBIC
#   method x criteria x nrep
# opt_criterias_lambdas: chosen lambda of each iteration, without cleaning up the error lambdas, according to AIC/BIC/EBIC
#   method x criteria x nrep
######################################################################################

filename = sprintf("results/datasets_N%0.0f_M%0.0f_nreps%0.0f.RData", N, M, nrep)
save.image(file = filename)

### Generate images

#source("results_gen_Coefs.R")
#source("results_gen_VS.R")

