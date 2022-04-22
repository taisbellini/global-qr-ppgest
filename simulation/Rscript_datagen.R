#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Inform parameters: N, nrep", call.=FALSE)
} else if (length(args)>0) {
  print(args)
}

N = as.numeric(args[1])
nrep = as.numeric(args[2])

#load libs

# Load aux files
source("core/data_generator.R")
source("core/utils.R")

# Monte Carlo

set.seed(205651)
start = Sys.time()

D= 30
X_array = array(0, dim=c(N, D, nrep))
Y_array = array(0, dim = c(N,nrep))
U_array = array(0, dim = c(N, nrep))

for (i in 1:nrep){

  # data gen function from data_generator.R
  data = data_gen(N=N)
  
  X_array[,,i] = data$X
  Y_array[,i] = data$y
  U_array[,i] = data$U

}

end = Sys.time()
time = end - start
print(time)

filename = sprintf("data_N%0.0f_nrep%0.0f.RData", N, nrep)
save(X_array, Y_array, U_array, file = filename)





