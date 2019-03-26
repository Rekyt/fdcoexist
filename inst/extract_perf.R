# Packages
library("purrr")
library("future.apply")
devtools::load_all()

# Given task
cli_args = commandArgs(TRUE)
job_task_id = as.numeric(cli_args[1])

number_of_sets_per_task = 1010


# Parameters used
list_A = c(0, 10^-(seq(1, 8, length.out = 6)))
list_k = seq(1, 1.5, length.out = 6)
list_B = list_A
list_H = seq(0, 1, length.out = 6)
n_seed = 30
n_patches = 25
n_species = 100
n_gen = 50
n_traits = 2

param_sets = list(
    run_n = seq(n_seed),
    A     = list_A,
    k     = list_k,
    B     = list_B,
    H     = list_H) %>%
    # Make all combinations but exclude cases where A > B
    cross(.filter = ~..2 > ..4)


# Return split sequence for a giving number of
f = function(a, b) {
    seq((a - 1) * b + 1, a * b, by = 1)
}

param_used = f(job_task_id, number_of_sets_per_task)

if (max(param_used) > length(param_sets)) {
    param_used = seq(min(param_used), length(param_sets), by = 1)
}

# Load data
print("Load data")
tictoc::tic()
var_param = readRDS(paste0("inst/job_data/var_param_bigmem_", min(param_used),
                           "_", max(param_used), "_data.Rds"))
tictoc::toc()

used_trait_list = readRDS("inst/job_data/traits_bigmem_30.Rds")

# Extract performances
plan(multicore, workers = 8)

tictoc::tic()
var_param_perfs = future_lapply(unlist(var_param, recursive = FALSE),
                          function(x) {

                              suppressMessages({devtools::load_all})

                              extract_performances_from_simul(x, used_trait_list,
                                                           TRUE)
                          })
tictoc::toc()

saveRDS(var_param_perfs, paste0("inst/job_data/perf_bigmem_", min(param_used),
                                "_", max(param_used),".Rds", compress = TRUE))
