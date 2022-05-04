# global-qr-ppgest
This repository contains the code used for the Master's Degree dissertation in PPGEst - UFRGS. Global Variable Selection for Quantile Regression.

## How to use

### Data Generation
#### Generate new data

Run the script passing the parameters for sample size and number of simulations.

```
Rscript Rscript_data_gen.R <N> <nrep>
```

It will generate a file called `data_N<>_nrep<>.RData` that you will use in the monte carlo. You can update the `set.seed()` line for different outcomes.

#### Get the simulation data used in thesis
Free github do not support very large files, so the simulated data X and Y are stored in: https://drive.google.com/file/d/1IWioTEkP5I3ujHJn6dYXXosAlJmtnl1L/view?usp=sharing under the name `data_N1000_nrep10000_official.RData`. Since it is a very big file, it occupies a lot of memory, so there is a subset of the first 1000 replications in file `data_N1000_nrep1000_official.RData`.

### Monte Carlo
Run the script passing the parameters for sample size, quantile levels grid size, batch of replications (start and end value), number of lambdas to be tested. The script will take the replications from `nrep_start` to `nrep_end` and run all the 6 methods evaluated in the thesis with the `N`, `M` and `lambdas` parameters. It will generate an file called `results_N<>f_M<>_from<>_to<>.RData` in the results/RData_outputs folder.

```
Rscript Rscript_monte_carlo.R <N> <M> <nrep_start> <nrep_end> <lambdas>
```



