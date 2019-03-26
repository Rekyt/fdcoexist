# Packages
library("purrr")
library("future.apply")
devtools::load_all()

# Parameters used
n_seed = 30
n_patches = 25
n_species = 100
n_traits = 2

# Get trait list
used_trait_list = map(seq(n_seed), function(given_seed) {
    set.seed(given_seed)
    given_traits = generate_cor_traits(n_patches, n_species, n_traits - 1,
                                       cor_coef = 0)
    set.seed(given_seed)
    cor_trait = given_traits = generate_cor_traits(n_patches, n_species,
                                                   n_traits - 1,
                                                   cor_coef = 0.3)
    set.seed(given_seed)
    negcor_trait = given_traits = generate_cor_traits(n_patches, n_species,
                                                      n_traits - 1,
                                                      cor_coef = -0.3)

    list(uncor  = given_traits,
         poscor = cor_trait,
         negcor = negcor_trait)
}) %>%
    set_names(nm = seq(n_seed))

saveRDS(used_trait_list, paste0("inst/job_data/traits_bigmem_", n_seed, ".Rds"),
        compress = TRUE)
