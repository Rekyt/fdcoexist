# Script to extract CWMs from simulations
# Author: Matthias Grenié <matthias dot grenie at gmail dot com>
# Packages ---------------------------------------------------------------------
library("dplyr")

# Functions --------------------------------------------------------------------
# Sequence from start to a given number even if last interval is smaller than
# expected
seqlast = function(from, to, by) {
    vec = do.call(what = seq, args = list(from, to, by))
    if (tail(vec, 1) != to) {
        return(c(vec, to))
    } else {
        return(vec)
    }
}

# Load initial data ------------------------------------------------------------
list_perf_files = list.files("inst/job_data/perfs_52920/", "all_perf_df_*",
                             full.names = TRUE, recursive = FALSE)

trait_df = readRDS("inst/job_data/perfs_52920/bigmem_trait_df.Rds") %>%
    mutate(seed = as.numeric(seed))

# Compute CWMs -----------------------------------------------------------------
lapply(list_perf_files[10], function(single_perf_file) {
    # Get the numbers in the perf file
    param_groups = basename(single_perf_file)

    number = regmatches(param_groups, gregexpr("[0-9]+_[0-9]+",
                                               param_groups))[[1]]

    cat("number: ", number, "\n")

    # Load performance data
    perf_data = readRDS(single_perf_file)

    # Get all combinations of paramters tested
    group_df = perf_data %>%
        distinct(h_fun, di_thresh, k, A, B, H, seed)

    # Split in groups of 1000
    groups = seqlast(1, nrow(group_df), 500)

    # Compute CWM over each group
    purrr::imap(groups[-length(groups)], function(x, y) {
        perf_data %>%
            semi_join(group_df %>%
                          slice(groups[y], groups[y + 1])) %>%
            inner_join(trait_df, by = c("seed", "species", "trait_cor")) %>%
            group_by(h_fun, di_thresh, k, A, B, H, R_scenar, A_scenar,
                     H_scenar, trait_cor, seed, patch) %>%
            summarise_at(vars(starts_with("trait")),
                         list(cwm = ~weighted.mean(., N150))) %>%
            saveRDS(paste0("inst/job_data/perfs_52920/cwm/cwm_groups_",
                           number, "_", groups[y], "_", groups[y + 1], ".Rds"),
                    compress = TRUE)
    })
})
