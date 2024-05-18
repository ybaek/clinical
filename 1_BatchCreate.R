###Start with a clean space
rm(list = ls())

###Server
server <- "hardac" # options are "hardac" or "dash"
if (server == "hardac") file_path <- "/data/mukherjeelab/Berchuck/final/paper1a/"
if (server == "dash") file_path <- "/data/berchucklab/berchuck/simulations/final/paper1a/"

###Set seed
set.seed(54)

###Batch settings
n_data <- 100
minibatch <- c(1, 5, 10)
delta <- seq(0, 1, 0.1)
sample_size <- c(10000, 100000)
settings <- expand.grid(sample_size, delta, minibatch)
settings$epsilon <- settings$Var3 / settings$Var1^(1 + settings$Var2)
settings$bound <- 1 / settings$Var1
settings <- settings[settings$epsilon < settings$bound, ]
rownames(settings) <- NULL
settings$epsilon <- settings$bound <- NULL
n_settings_unique <- nrow(settings)
settings <- data.frame(setting = 1:n_settings_unique, settings)
colnames(settings)[2:4] <- c("n", "delta", "minibatch")
settings <- data.frame(settings[rep(1:n_settings_unique, each = n_data), ], data = rep(1:n_data, n_settings_unique))
rownames(settings) <- NULL
n_settings <- nrow(settings)
settings$epsilon <- settings$minibatch / settings$n^(1 + settings$delta)
settings$nsims <- 100 / settings$epsilon

###Printing functions
dprint <- function(x) cat(paste0(x, "\n"), file = file_out, append = TRUE)
skip <- function() cat("\n", file = file_out, append = TRUE)

###Loop over settings
for (s in 1:n_settings) {

  ###Simulation settings
  setting <- settings$setting[s]
  n <- settings$n[s]
  data <- settings$data[s]
  delta <- settings$delta[s]
  S <- settings$minibatch[s]

  ###Create output file
  file_out <- paste0(file_path, "scripts/", s, ".R")
  cat(file = file_out)
  if (server == "hardac") {
    dprint("###Set the R library path")
    dprint(".libPaths(c(.libPaths(), \"/home/sib2/R\"))")
    skip()
  }
  dprint("###Load packages")
  dprint("library(sglmm)")
  skip()
  dprint("###Set seed")
  dprint(paste0("seed <- ", data))
  dprint("set.seed(seed)")
    skip()
  dprint("###Run simulation")
  dprint(paste0("out <- paper_sim1_Rcpp(S = ", S, ", delta = ", delta, ", n = ", n, ", t_real = 1)"))
  dprint(paste0("source(file = \"", file_path, "scripts/source_stan.R\")"))
  dprint(paste0("save(out, file = \"", file_path, "output/", s, ".RData\")"))
    skip()
  dprint("###End session")
  dprint("q(save = \"no\")")
    skip()
  
###End loop over settings
}

###Run everything in an array!
file_out <- paste0(file_path, "scripts/RunAll.sh")
cat(file = file_out)
dprint("#!/bin/bash")
skip()
dprint("#SBATCH -N 1")
dprint("#SBATCH -t 10-00:00:00 # Wall time limit in DD-HH:MM:SS")
dprint("#SBATCH -J paper1")
dprint("#SBATCH -o paper1.out")
dprint("#SBATCH -e paper1.err")
if (server == "dash") dprint("#SBATCH --account=berchucklab")
dprint("#SBATCH --mem=4G")
dprint(paste0("#SBATCH --array=1-", n_settings, "%", 100))
dprint("#SBATCH --mail-type=FAIL")
dprint("#SBATCH --mail-user=sib2@duke.edu")
skip()
if (server == "hardac") {
  dprint("module load Anaconda3/2019.10-gcb02")
  # dprint("conda activate /data/common/shared_conda_envs/miniconda3/envs/R-4.0.2")
  dprint("conda activate R")
}
if (server == "dash") {
  dprint("source /sched/anaconda3/etc/profile.d/conda.sh")
  dprint("conda activate R")
}
dprint(paste0("srun R CMD BATCH \"", file_path, "scripts/$SLURM_ARRAY_TASK_ID.R\""))
dprint("conda deactivate")