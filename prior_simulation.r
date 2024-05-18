source("prior_model.r")
library(dplyr)
library(ggplot2)
library(tidyr)

set.seed(20435619)

N <- 10000
Pr_obs <- 0.1
N_obs <- round(N * Pr_obs)
N_mis <- N - N_obs
P <- 100
total_count <- 40
thres_grid <- seq(0, total_count - 1)
beta <- rnorm(P) / sqrt(P)
gamma <- sqrt(.8) * beta + rnorm(P, 0, sqrt(.2 / P))
beta0_grid <- seq(-1, 3)
gamma0_grid <- seq(-2, 2)
beta_R_grid <- seq(.1, 1, by = .1)
design <- expand.grid(
    betaR = beta_R_grid,
    gamma0 = gamma0_grid,
    beta0 = beta0_grid
)
##################################################

#   Monte Carlo to address conditioning issues   #

##################################################
design_base <- design
pnv_dfs <- list()
for (i in 1:dim(design_base)[1]) {
    pars <- as.numeric(design_base[i, ])
    beta_R <- pars[1]
    gamma0 <- pars[2]
    beta0 <- pars[3]
    X <- matrix(rnorm(N * P), N, P) / sqrt(P)
    R <- generate_R(X, gamma, gamma0)
    Y <- generate_Y(X, R, total_count, gamma, gamma0, beta, beta0, beta_R)
    ppvs <- ppvs_stat <- numeric(total_count)
    npvs <- npvs_stat <- numeric(total_count)
    sens_stat <- spec_stat <- numeric(total_count)
    senss <- specs <- numeric(total_count)
    for (j in 1:total_count) {
        thres <- thres_grid[j]
        Z <- Y > thres
        ppvs[j] <- sum(Z == 1 & R == 1) / sum(R == 1)
        npvs[j] <- sum(Z == 0 & R == 0) / sum(R == 0)
        senss[j] <- sum(Z == 1 & R == 1) / sum(Z == 1)
        specs[j] <- sum(Z == 0 & R == 0) / sum(Z == 0)
        # Compare with other approximations
        approx <- generate_stats(
            X, total_count, beta, beta0, beta_R, gamma, gamma0, thres
        )
        ppvs_stat[j] <- approx[2]
        npvs_stat[j] <- approx[3]
        sens_stat[j] <- approx[4]
        spec_stat[j] <- approx[5]
    }
    pnv_dfs[[i]] <- data.frame(
        ppv_mc1 = ppvs,
        npv_mc1 = npvs,
        sens_mc1 = senss,
        spec_mc1 = specs,
        ppv_mc2 = ppvs_stat,
        npv_mc2 = npvs_stat,
        sens_mc2 = sens_stat,
        spec_mc2 = spec_stat,
        thres = thres_grid,
        betaR = rep(beta_R, total_count),
        gamma0 = rep(gamma0, total_count),
        beta0 = rep(beta0, total_count)
    )
}
pnvs <- dplyr::bind_rows(pnv_dfs)
#
# 1. How do PPV and NPV look?
gpn <- ggplot(pnvs, aes(x = ppv_mc1, y = npv_mc1)) +
    geom_line(aes(color = betaR, group = betaR), alpha = .5) +
    facet_grid(beta0 ~ gamma0, labeller = label_both) +
    labs(x = "PPV", y = "NPV")
ggsave("figures/ppv-npv-dense.jpeg",
    plot = gpn,
    device = "jpeg",
    width = 20,
    height = 20,
    units = "cm"
)
#
# 2. How do we do on Sens/Spec in terms of Monte Carlo?
gss1 <- pnvs %>%
    ggplot(aes(x = 1 - spec_mc1, y = sens_mc1)) +
    geom_line(aes(color = betaR, group = betaR), alpha = .5) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_grid(beta0 ~ gamma0, labeller = label_both) +
    labs(x = "1 - Spec", y = "Sens", title = "Fully sampled")
gss2 <- pnvs %>%
    ggplot(aes(x = 1 - spec_mc2, y = sens_mc2)) +
    geom_line(aes(color = betaR, group = betaR), alpha = .5) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_grid(beta0 ~ gamma0, labeller = label_both) +
    labs(x = "1 - Spec", y = "Sens", title = "Only X sampled")
ggsave("figures/sens-spec-dense1.jpeg",
    plot = gss1,
    device = "jpeg",
    width = 20,
    height = 20,
    units = "cm"
)
ggsave("figures/sens-spec-dense2.jpeg",
    plot = gss2,
    device = "jpeg",
    width = 20,
    height = 20,
    units = "cm"
)
# !
# We get some NaNs.
# !
nan_tbl1 <- pnvs %>%
    mutate(
        nan_sens = is.nan(sens_mc1), nan_spec = is.nan(spec_mc1)
    ) %>%
    group_by(betaR, beta0, gamma0, thres) %>%
    summarise(nan_sens = sum(nan_sens), nan_spec = sum(nan_spec)) %>%
    pivot_longer(
        cols = c("nan_sens", "nan_spec"),
        names_to = "Stat", values_to = "NaNs"
    )
nan_tbl2 <- pnvs %>%
    mutate(
        nan_sens = is.nan(sens_mc2), nan_spec = is.nan(spec_mc2)
    ) %>%
    group_by(betaR, beta0, gamma0, thres) %>%
    summarise(nan_sens = sum(nan_sens), nan_spec = sum(nan_spec)) %>%
    pivot_longer(
        cols = c("nan_sens", "nan_spec"),
        names_to = "Stat", values_to = "NaNs"
    )
gnan1 <- nan_tbl1 %>%
    ggplot(aes(x = thres, y = NaNs)) +
    geom_bar(aes(group = Stat, fill = Stat), stat = "identity") +
    facet_grid(beta0 ~ gamma0, labeller = label_both) +
    labs(fill = "Sens/Spec", title = "Fully Sampled")
gnan2 <- nan_tbl2 %>%
    ggplot(aes(x = thres, y = NaNs)) +
    geom_bar(aes(group = Stat, fill = Stat), stat = "identity") +
    facet_grid(beta0 ~ gamma0, labeller = label_both) +
    labs(fill = "Sens/Spec", title = "Only X Sampled")
ggsave("figures/nans-dense1.jpeg",
    plot = gnan1,
    device = "jpeg",
    width = 20,
    height = 20,
    units = "cm"
)
ggsave("figures/nans-dense2.jpeg",
    plot = gnan2,
    device = "jpeg",
    width = 20,
    height = 20,
    units = "cm"
)
###################################################
#                                                 #
#        To be figured out                        #
#                                                 #
###################################################
# 2. Sample missingness mechanism
alpha <- sqrt(.8) * beta + rnorm(P, 0, sqrt(.2 / P))
alpha0 <- 0
alpha_YR <- 0
alpha_Y_grid <- seq(-1, 1, by = .5) / total_count * 10
alpha_R_grid <- seq(-1, 1, by = .5)
design_mis <- expand.grid(
    alpha_Y = alpha_Y_grid,
    alpha_R = alpha_R_grid
)
o_dfs <- list()
for (i in 1:dim(design)[1]) {
    pars <- as.numeric(design[i, ])
    beta_R <- pars[1]
    gamma0 <- pars[2]
    beta0 <- pars[3]
    X <- matrix(rnorm(N * P), N, P) / sqrt(P)
    R <- generate_R(X, gamma, gamma0)
    Y <- generate_Y(X, R, total_count, gamma, gamma0, beta, beta0, beta_R)
    for (j in 1:dim(design_mis)[1]) {
        pars <- as.numeric(design_mis[j, ])
        alpha_Y <- pars[1]
        alpha_R <- pars[2]
        O <- generate_O(X, Y, R, alpha, alpha0, alpha_Y, alpha_R, alpha_YR)
        o_dfs[[(i - 1) * dim(design_mis)[1] + j]] <- data.frame(
            no_mis = mean(R == 0 & O == 0),
            yes_mis = mean(R == 1 & O == 0),
            no_obs = mean(R == 0 & O == 1),
            yes_obs = mean(R == 1 & O == 1),
            alpha_Y = alpha_Y,
            alpha_R = alpha_R,
            beta0 = beta0,
            beta_R = beta_R,
            gamma0 = gamma0
        )
    }
}
o_dfs <- dplyr::bind_rows(o_dfs)
o_dfs <- o_dfs %>%
    pivot_longer(
        cols = c("no_mis", "yes_mis", "no_obs", "yes_obs"),
        names_to = "Stat", values_to = "Props"
    )
gmis1 <- o_dfs %>%
    filter(beta_R < 0.55 & beta_R > 0.45 & Stat == "yes_mis") %>%
    ggplot(aes(x = beta0, y = gamma0)) +
    geom_tile(aes(fill = Props)) +
    facet_grid(alpha_Y ~ alpha_R, labeller = label_both) +
    labs(title = "Missing with R = 1 (beta_R = 0.5)")
gmis2 <- o_dfs %>%
    filter(beta_R < 0.55 & beta_R > 0.45 & Stat == "no_mis") %>%
    ggplot(aes(x = beta0, y = gamma0)) +
    geom_tile(aes(fill = Props)) +
    facet_grid(alpha_Y ~ alpha_R, labeller = label_both) +
    labs(title = "Missing with R = 0 (beta_R = 0.5)")
ggsave("figures/missing_test1.jpeg",
    plot = gmis1,
    device = "jpeg",
    width = 20,
    height = 20,
    units = "cm"
)
ggsave("figures/missing_test2.jpeg",
    plot = gmis2,
    device = "jpeg",
    width = 20,
    height = 20,
    units = "cm"
)
