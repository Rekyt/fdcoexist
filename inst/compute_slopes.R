# Script to get CWM slopes
# Packages ---------------------------------------------------------------------
library("dplyr")
library("future.apply")

# Functions --------------------------------------------------------------------
seqlast = function(from, to, by) {
    vec = do.call(what = seq, args = list(from, to, by))
    if (tail(vec, 1) != to) {
        return(c(vec, to))
    } else {
        return(vec)
    }
}

# Load Data --------------------------------------------------------------------
all_cwm = readRDS("inst/job_data/perf_2fd398/all_cwm.Rds")

all_groups = all_cwm %>%
    distinct(k, A, B, H, h_fun, di_thresh, R_scenar, A_scenar, H_scenar,
             trait_cor, seed)

n_groups = nrow(all_groups)

group_list = seqlast(1, n_groups, 1e5)



# Cut CWM file into slices -----------------------------------------------------
lapply(seq_along(group_list[-length(group_list)]),
       function(group_index) {
           all_cwm %>%
               semi_join(all_groups %>%
                             slice(group_list[group_index],
                                   group_list[group_index + 1]),
                         by = c("k", "A", "B", "H", "h_fun", "di_thresh",
                                "R_scenar", "A_scenar", "H_scenar",
                                "trait_cor", "seed")) %>%
               saveRDS(paste0("inst/job_data/perf_2fd398/cwm_slice_",
                                       group_list[group_index], "_",
                                       group_list[group_index + 1], ".Rds"),
                       compress = TRUE)
})

# Compute CWM-Env linear models on slices  -------------------------------------
future::plan(future::multiprocess, workers = 62)

cwm_slices = list.files("inst/job_data/perf_2fd398", "cwm_slice_*")

future_lapply(
    cwm_slices,
    function(cwm_slice) {

        slice_group = regmatches(cwm_slice, regexpr("[0-9]+_[0-9]+", cwm_slice))

        cwm_mod = cwm_slice %>%
            select(-trait2_cwm, -trait1_cwv, -trait2_cwv, -species_rich) %>%
            mutate(trait1_cwm = ifelse(is.na(trait1_cwm), 0, trait1_cwm)) %>%
            tidyr::nest(trait1_cwm, patch) %>%
            mutate(lm_mod = purrr::map(data, ~lm(trait1_cwm ~ patch,
                                                 data = .x)),
                   lm_sum = purrr::map(lm_mod, broom::tidy))

        saveRDS(cwm_mod, paste0("inst/job_data/perf_2fd398/cwm_mod_",
                                slice_group, ".Rds"),
                compress = TRUE)

        cwm_mod %>%
            unnest(lm_sum) %>%
            filter(estimate == "patch") %>%
            saveRDS(paste0("inst/job_data/perf_2fd398/cwm_slopes_",
                           slice_group, ".Rds"),
                    compress = TRUE)
    }
)
