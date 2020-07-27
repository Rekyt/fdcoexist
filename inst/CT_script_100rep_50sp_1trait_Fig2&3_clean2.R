# FDcoexist June 4 2029
# C Tucker
# Generate species level responses to site conditions. Only moderate value of A,
# H or both used.
# replicates of each required - currently ~100
# For plots 2 & 3 (update with st err)

# Load Packages ----------------------------------------------------------------
library(dplyr)
library(cowplot)
library(tidyr)
library(RColorBrewer)
library(sads)
devtools::load_all()  # Load all functions from package in local repo



# Initial Parameters -----------------------------------------------------------
set.seed(1)      # set random seed
n_patches <- 25  # number of patches
n_species <- 50  # number of species
n_gen <- 50      # number of time steps
n_traits <- 2    # number of traits
init_pop <- 50   # number of individuals at t=0 for each species
d <- 0.05        # dispersal parameter
width <- 5       # standard deviation of the Gaussian environmental filtering
K <- 100         # carrying capacity (per species per patch)

# Initial population matrix, 50 individuals of each species
composition <- array(NA, dim = c(n_patches, n_species, n_gen*20),
                     dimnames = list(paste0("patches", seq(n_patches)),
                                     paste0("species", seq(n_species)),
                                     paste0("time", seq(n_gen*20))))

composition[, , 1] <- init_pop


# Trait scenario with one trait contributing to all processes
trait_scenar <- data.frame(trait            = "trait1",
                           growth_weight    = 1,
                           compet_weight    = 1,
                           hierarchy_weight = 1)

# One trait, all four scenarios with X replicates and random trait var.
# Trait value (uncorrelated traits)
uncor_traits <- generate_cor_traits(n_patches, n_species, n_traits - 1,
                                    cor_coef = 0)
uncor_traits <- uncor_traits[, 1, drop = FALSE]
trait_comb <- seq(1:100)

# Generate truly random trait distributions (trait 1 varies)
for (i in trait_comb[-1]) {
    uncor_traits2 <- generate_cor_traits_rand(n_patches, n_species,
                                              n_traits - 1, cor_coef = 0)
    uncor_traits2 <- uncor_traits2[, 1, drop = FALSE]
    uncor_traits <- cbind(uncor_traits, uncor_traits2)
}

# One trait contributing to all the processes
# Parameter values - reduced for speed
list_k <- c(2)
list_A <- c(0, 1e-3)
list_H <- c(0, 1e-3)

# data.frame with all combinations of parameter values
comb <- data.frame(expand.grid(list_k, list_A, list_H, trait_comb))
colnames(comb) <- c("k", "A", "H", "trait_comb")
comb <- distinct(comb) # remove duplicates


# Computation & Outputs --------------------------------------------------------
eq <- c()
mis_i <- list()
tra_env <- list()

# growth rates only estimated for final time step
j = 2*n_gen

for (i in 1:nrow(comb)) {

    # Extract simulation trait data.frame
    traits <- data.frame(trait1 = uncor_traits[,comb[i, "trait_comb"]])
    simul_i <- multigen(
        traits = traits, trait_weights = trait_scenar, env = 1:n_patches,
        time = 2*n_gen, species = n_species, patches = 25,
        composition = composition,
        A = comb[i, "A"], B = 1e-7, d = d, k = comb[i, "k"],
        H = comb[i, "H"],
        width = rep(width, n_patches), h_fun = "+", di_thresh = 24, K = 100)


    ## Data frame of last time step with traits
    tra_env_j <- reshape2::melt(simul_i$compo[, , j])
    colnames(tra_env_j) <- c("site", "sp", "ab")

    sp_tra <- data.frame(sp = rownames(uncor_traits),
                         tra = uncor_traits[, comb[i, "trait_comb"]])
    tra_env_j <- suppressWarnings(left_join(tra_env_j, sp_tra, by = "sp"))

    ## Compute CWM
    tra_env_j <- tra_env_j %>%
        group_by(site) %>%
        mutate(time = j,
               k = comb[i,"k"], A = comb[i, "A"], H = comb[i, "H"],
               A_H = comb[i, "trait_comb"],
               cwm = weighted.mean(tra, w = ab, na.rm=TRUE)) %>%
        as.data.frame() %>%
        mutate(env = as.integer(gsub("patches", "", x = site)))

    ## Compute growth rates
    max_time <- n_gen
    gw <- extract_growth_rates(simul=simul_i, chosen_time = j)

    # Merge growth rates with the CWM data frame
    tra_env_j <- left_join(tra_env_j, gw, by = c("env" = "patch", "sp" = "species"))

    # Add trait weighted by the different growth rates per site
    # (scaling before)
    tra_env_j <- tra_env_j %>%
        group_by(site) %>% # filter(site == "patches15") %>%
        mutate(avg_growth_rate = scales::rescale(avg_growth_rate,
                                                 c(0, 1)),
               max_growth_rate = scales::rescale(max_growth_rate,
                                                 c(0, 1)),
               int_growth_rate = scales::rescale(int_growth_rate,
                                                 c(0, 1))) %>%
        mutate(w_avg_gw = wtd_mean(tra, w = avg_growth_rate, na.rm=TRUE),
               w_max_gw = wtd_mean(tra, w = max_growth_rate, na.rm=TRUE),
               w_int_gw = wtd_mean(tra, w = int_growth_rate, na.rm=TRUE)) %>%
        distinct(site, .keep_all = TRUE) %>%
        select(-sp, -final_abundance, -avg_growth_rate,
               -max_growth_rate, -int_growth_rate) %>%
        as.data.frame()

    ## Bind patch level results
    tra_env[[i]] <- tra_env_j

    ## Species mismatches
    all_growth <- r_env_CT(simul_i, sp1 = n_species, n_patches = n_patches,
                           time = 50)

    ##species level mismatches
    matches <- matrix(NA, ncol=17, nrow=n_species)
    for (z in 1:n_species){
        allout <- extract_mismatchesCT(x=all_growth, z=z)

        sub <- na.omit(gw[which(gw$sp == paste("species",z, sep="")), ])
        if(nrow(sub)>0){
            finalabund <- sub[which.max(sub$final_abundance), "patch"]
            avg_growth_rate <- sub[which.max(sub$avg_growth_rate), "patch"]
            max_growth_rate <- sub[which.max(sub$max_growth_rate), "patch"]
            int_growth_rate <- sub[which.max(sub$int_growth_rate), "patch"]
            finalabund1 <- max(sub $final_abundance, na.rm=TRUE)
            avg_growth_rate1 <- max(sub$avg_growth_rate, na.rm=TRUE)
            max_growth_rate1 <- max(sub$max_growth_rate, na.rm=TRUE)
            int_growth_rate1 <- max(sub$int_growth_rate, na.rm=TRUE)

            matches[z,] <- c(z, allout, finalabund, finalabund1, avg_growth_rate, avg_growth_rate1, max_growth_rate, max_growth_rate1, int_growth_rate, int_growth_rate1)
        }else{

            matches[z,] <- c(z, allout, NA, NA, NA, NA, NA, NA, NA, NA)
        }
    }

    colnames(matches) <- c("Species", "TrueRPatch", "TrueR", "ObsRPatch",
                           "ObsR", "TrueAbPatch", "TrueAb", "ObsAbPatch",
                           "ObsAb", "finalabundPatch", "finalabund",
                           "avgGRPatch", "avgGR", "maxGRPatch", "maxGR",
                           "intGRPatch", "intGR")

    matches <- cbind(matches, comb[rep(i, nrow(matches)), ], traits)

    mis_i[[i]] <- matches

    cat(100*i/nrow(comb), "\t")
}

#Save files
#tra_env100 <- tra_env
#mis_i100 <- mis_i


# Mismatches by Species --------------------------------------------------------

allmis <- bind_rows(mis_i100)

# Calculate all the various mismatches. Also super messy.

# Final abundance is true env only, vs obs abundance (MisAbP)
# prop change in obs ab from TRUE
allmis$MisAbP <- 100*(allmis$finalabund - allmis$ObsAb)/n_patches
#by patch match
allmis$MisAbPatchP <- 100*(allmis$finalabundPatch - allmis$ObsAbPatch)/n_patches

#Obs inst GR as measure
# prop change in obs instR patch
allmis$MisinstRP <- 100*(allmis$intGR - allmis$ObsR)/n_patches
# patch match
allmis$MisinstRPatchP <- 100*(allmis$intGRPatch - allmis$ObsRPatch)/n_patches

# Avg GR as measure
allmis$MisavgGRP <- 100*(allmis$avgGR - allmis$ObsR)/n_patches #prop change in obs instR
allmis$MisavgGRPatchP <- 100*(allmis$avgGRPatch - allmis$ObsRPatch)/n_patches

# Max GT as a measure
allmis$MismaxGRP <- 100*(allmis$maxGR - allmis$ObsR)/n_patches
allmis$MismaxGRPatchP <- 100*(allmis$maxGRPatch - allmis$ObsRPatch)/n_patches

allmis$comb <- paste(allmis$k, allmis$A, allmis$H, sep="_")

 #ignore warnings
allmisA <- aggregate(allmis, by=list(allmis$Species, allmis$comb), "mean", na.rm=TRUE)
allmisA$comb <- allmisA$Group.2

unique(allmisA$comb)
#subset into distinct scenarios
nocomp <- allmisA[which(allmisA$comb %in% c("2_0_0")), ]
Acomp <- allmisA[which(allmisA$comb %in% c("2_0.001_0")), ]
Hcomp <- allmisA[which(allmisA$comb %in% c("2_0_0.001")), ]
allcomp <- allmisA[which(allmisA$comb %in% c("2_0.001_0.001")), ]


# Plots of Mismatches by Species -----------------------------------------------
##mismatch plots, by species
par(mfrow=c(1,2))
plot(nocomp$Species~ nocomp$MisAbPatchP, cex=0.3, xlim=c(-25,25),
     xlab="Mismatch from true topt (% of gradient)", ylab="Species",
     main="+Limiting similarity", col="grey", type="l")
abline(h=c(1:50), col="grey", lwd=0.2)
points(Acomp$Group.1~Acomp$MisAbPatchP, cex=0.75, col="red")
points((Acomp$Group.1) ~ (Acomp$MisinstRPatchP), cex=0.75, pch=16, col="orange")
arrows(x0= Acomp$MisAbPatchP, x1=Acomp$MisinstRPatchP, y0=seq(1:50), y1=seq(1:50),code=0)
points((Acomp$Species) ~ (Acomp$MismaxGRPatchP), cex=0.75, pch=09, col="blue")
points((Acomp$Species) ~ (Acomp$MisavgGRPatchP), cex=0.75, pch=4, col="purple")
arrows(x0= Acomp$MismaxGRPatchP, x1=Acomp$MisavgGRPatchP, y0=seq(1:50), y1=seq(1:50),code=0)
arrows(x0= Acomp$MisPatchAbP, x1=Acomp$MisavgGRPatchP, y0=seq(1:50), y1=seq(1:50),code=0)

points(y=c(50), x=10,  col="red")
text(y=c(50), x=20, cex=0.75,"Obs abundance")
points(y=c(48), x=10, pch=16, col="orange")
text(y=c(48), x=20, cex=0.75,"Inst. growth rate")
points(y=c(46), x=10,  pch=9, col="blue")
text(y=c(46), x=20, cex=0.75,"Max growth rate")
points(y=c(44), x=10,  pch=4, col="purple")
text(y=c(44), x=20, cex=0.75,"Avg growth rate")
#

#H competition
plot(nocomp$Species~ nocomp$MisAbPatchP, cex=0.3, xlim=c(-25,25), xlab="Mismatch from true topt (% of gradient)", ylab="Species", main="+Hierarchical competition", type="l", col="grey")
#points(I(nocomp$Species) ~ (nocomp$MisPatchinstRP), cex=0.3, col="grey")
#abline(v=c(0), col="grey", lwd=1)
abline(h=c(1:50), col="grey", lwd=0.2)
points(Hcomp$Group.1~Hcomp$MisAbPatchP, cex=0.75, col="red")
points((Hcomp$Group.1) ~ (Hcomp$MisinstRPatchP), cex=0.75, pch=16, col="orange")
arrows(x0= Hcomp$MisAbPatchP, x1=Hcomp$MisinstRPatchP, y0=seq(1:50), y1=seq(1:50),code=0)
points((Hcomp$Species) ~ (Hcomp$MismaxGRPatchP), cex=0.75, pch=09, col="blue")
points((Hcomp$Species) ~ (Hcomp$MisavgGRPatchP), cex=0.75, pch=4, col="purple")
arrows(x0= Hcomp$MismaxGRPatchP, x1=Hcomp$MisavgGRPatchP, y0=seq(1:50), y1=seq(1:50),code=0)
arrows(x0= Hcomp$MisAbPatchP, x1=Hcomp$MisavgGRPatchP, y0=seq(1:50), y1=seq(1:50),code=0)

#all comp
#plot(nocomp$Species~ nocomp$MisPatchAbP, cex=0.3, xlim=c(-25,25), xlab="Mismatch from true topt (% of gradient)", ylab="Species", main="All processes", col="grey", type="l")
#abline(h=c(1:50), col="grey", lwd=0.2)
#points(allcomp$Group.1~allcomp$MisPatchAbP, cex=0.75, col="red")
#points((allcomp$Group.1) ~ (allcomp$MisPatchinstRP), cex=0.75, pch=16, col="orange")
#arrows(x0= allcomp$MisPatchAbP, x1=allcomp$MisPatchinstRP, y0=seq(1:50), y1=seq(1:50),code=0)
#points((allcomp$Species) ~ (allcomp$MismaxGRPatchP), cex=0.75, pch=09, col="blue")
#points((allcomp$Species) ~ (allcomp$MisavgGRPatchP), cex=0.75, pch=4, col="purple")
#arrows(x0= allcomp$MismaxGRPatchP, x1=allcomp$MisavgGRPatchP, y0=seq(1:50), y1=seq(1:50),code=0)
#arrows(x0= allcomp$MisPatchAbP, x1=allcomp$MisavgGRPatchP, y0=seq(1:50), y1=seq(1:50),code=0)


# Correlation of Mismatches ----------------------------------------------------
###tables of correlations
#allcomp
cor(allcomp[,c("MisAbPatchP", "MisinstRPatchP", "MismaxGRPatchP",
               "MisavgGRPatchP")],
    method="spearman")

# MisPatchAb MisPatchinstR MismaxGRPatch MisavgGRPatch
# MisPatchAb     1.0000000     0.4164496     0.6729627     0.8061727
# MisPatchinstR  0.4164496     1.0000000     0.3293066     0.5395412
# MismaxGRPatch  0.6729627     0.3293066     1.0000000     0.8102633
# MisavgGRPatch  0.8061727     0.5395412     0.8102633     1.0000000

#H comp
cor(Hcomp[,c("MisAbPatchP", "MisinstRPatchP", "MismaxGRPatchP",
             "MisavgGRPatchP")],
    method="spearman")
# MisPatchAb MisPatchinstR MismaxGRPatch MisavgGRPatch
# MisPatchAb     1.0000000    0.40114217    0.10229505    0.49460843
# MisPatchinstR  0.4011422    1.00000000    0.04985318    0.82697422
# MismaxGRPatch  0.1022950    0.04985318    1.00000000   -0.09408434
# MisavgGRPatch  0.4946084    0.82697422   -0.09408434    1.00000000

cor(Acomp[,c("MisAbPatchP", "MisinstRPatchP", "MismaxGRPatchP",
             "MisavgGRPatchP")], method="spearman")
# MisPatchAb MisPatchinstR MismaxGRPatch MisavgGRPatch
# MisPatchAb     1.0000000     0.7729351     0.6976272     0.9281220
# MisPatchinstR  0.7729351     1.0000000     0.4811728     0.7313187
# MismaxGRPatch  0.6976272     0.4811728     1.0000000     0.7848516
# MisavgGRPatch  0.9281220     0.7313187     0.7848516     1.0000000

#no comp
cor(nocomp[,c("MisAbPatchP", "MisinstRPatchP", "MismaxGRPatchP",
              "MisavgGRPatchP")], method="spearman")
# MisPatchAb MisPatchinstR MismaxGRPatch MisavgGRPatch
# MisPatchAb             1            NA            NA            NA
# MisPatchinstR         NA             1            NA            NA
# MismaxGRPatch         NA            NA             1            NA
# MisavgGRPatch         NA            NA            NA             1


# Mismatches by patch ----------------------------------------------------------
cwmdat <- bind_rows(tra_env100)
cwmdat$comb <- paste(cwmdat$k, cwmdat$A, cwmdat$H, sep="_")
unique(cwmdat$comb)
cwmdat$site2 <- as.numeric(gsub("[[:alpha:]]","", cwmdat$site))
#ignore all warning
##CMW by scenarios
nocompCWM <- subset(cwmdat, cwmdat$A==0 & cwmdat$H==0)
nocompCWM <- na.omit(nocompCWM)
nocompCWMA <- aggregate(nocompCWM, by=list(nocompCWM$site), mean, na.rm=TRUE)
nocompCWMV <- aggregate(nocompCWM[,9:15], by=list(nocompCWM$site), sd, na.rm=TRUE)
#
AcomponlyCWM <- subset(cwmdat, cwmdat$A>0 & cwmdat$H==0)
AcomponlyCWM <- na.omit(AcomponlyCWM)
AcomponlyCWMA <- aggregate(AcomponlyCWM, by=list(AcomponlyCWM$site), mean, na.rm=TRUE)
AcomponlyCWMV <- aggregate(AcomponlyCWM[,9:15], by=list(AcomponlyCWM$site), sd, na.rm=TRUE)
#
HcomponlyCWM <- subset(cwmdat, cwmdat$H>0 & cwmdat$A==0)
HcomponlyCWM <- na.omit(HcomponlyCWM)
HcomponlyCWMA <- aggregate(HcomponlyCWM, by=list(HcomponlyCWM$site), mean, na.rm=TRUE)
HcomponlyCWMV <- aggregate(HcomponlyCWM[,9:15], by=list(HcomponlyCWM$site), sd, na.rm=TRUE)

#
allcomponlyCWM <- subset(cwmdat, cwmdat$H>0 & cwmdat$A==0)
allcomponlyCWM <- na.omit(allcomponlyCWM)
allcomponlyCWMA <- aggregate(allcomponlyCWM, by=list(allcomponlyCWM$site), mean, na.rm=TRUE)
allcomponlyCWMV <- aggregate(allcomponlyCWM[,9:15], by=list(allcomponlyCWM$site), sd, na.rm=TRUE)

#Plot
par(mfrow=c(1,3), mar=c(4, 4, 2, 2))
plot(I((cwm-site2)/25)~site2, data= nocompCWMA, xlim=c(0,25),ylim=c(-0.15, 0.15), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="No competition (Env filtering only)", col="red")
abline(h=0, col="grey80")
points(I((w_int_gw-site2)/25)~site2, data= nocompCWMA, xlim=c(0,25), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="No competition (Env filtering only)", col="orange", pch=16)
points(I((w_max_gw-site2)/25)~site2, data= nocompCWMA, xlim=c(0,25), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="No competition (Env filtering only)", col="blue", pch=9)
points(I((w_avg_gw-site2)/25)~site2, data= nocompCWMA, xlim=c(0,25), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="No competition (Env filtering only)", col="purple", pch=4)
#
plot(I((cwm-site2)/25)~site2, data= AcomponlyCWMA, xlim=c(0,25),ylim=c(-0.15, 0.15), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="+ Limiting similarity", col="red")
abline(h=0, col="grey80")
points(I((w_int_gw-site2)/25)~site2, data= AcomponlyCWMA, xlim=c(0,25), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="No competition (Env filtering only)", col="orange", pch=16)
points(I((w_max_gw-site2)/25)~site2, data= AcomponlyCWMA, xlim=c(0,25), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="No competition (Env filtering only)", col="blue", pch=9)
points(I((w_avg_gw-site2)/25)~site2, data= AcomponlyCWMA, xlim=c(0,25), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="No competition (Env filtering only)", col="purple", pch=4)
#
plot(I((cwm-site2)/25)~site2, data= HcomponlyCWMA, xlim=c(0,25),ylim=c(-0.15, 0.15), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="+ Hierarchical competition", col="red")
abline(h=0, col="grey80")
points(I((w_int_gw-site2)/25)~site2, data= HcomponlyCWMA, xlim=c(0,25), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="No competition (Env filtering only)", col="orange", pch=16)
points(I((w_max_gw-site2)/25)~site2, data= HcomponlyCWMA, xlim=c(0,25), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="No competition (Env filtering only)", col="blue", pch=9)
points(I((w_avg_gw-site2)/25)~site2, data= HcomponlyCWMA, xlim=c(0,25), type="b", ylab="Obs CWM - expected (%)", xlab="Patch", main="No competition (Env filtering only)", col="purple", pch=4)

##These need to have error bars added!
save.image(file="/Users/tuckerca/Dropbox/Documents/Caroline/fdcoexist/2020/Data/100rep_50sp_1trait_scenarios.RData")
