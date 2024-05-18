rm(list = ls())
source("prior_model.r")
library(dplyr)
library(ggplot2)
library(tidyr)

set.seed(205619765)

N <- 10000
Pr_obs <- 0.1
N_obs <- round(N * Pr_obs)
N_mis <- N - N_obs
P <- 100
total_count <- 40
thres <- 14
beta <- rnorm(P) / sqrt(P)
alpha <- rnorm(P) / sqrt(P)
gamma <- sqrt(.8) * beta + rnorm(P, 0, sqrt(.2 / P))
design <- matrix(
    c(
        1.5, 1.5, -1.4, # gamma0
        1, -0.5, 0, # beta0
        0.2, 0.05, 0.01 # beta_R
    ),
    nrow = 3, ncol = 3
)
D <- dim(design)[1]
design_mis <- matrix(
    c(
        0.1, 0, # alpha0
        0.01, -0.04, # alpha_Y
        0, -1 # alpha_R
    ),
    nrow = 2, ncol = 3
)
Dd <- dim(design_mis)[1]
pr <- sens <- spec <- numeric(D)
pr_se <- sens_se <- spec_se <- numeric(D)
sens_obs <- spec_obs <- matrix(0, Dd, D)
sens_mis <- spec_mis <- matrix(0, Dd, D)
sens_obs_se <- spec_obs_se <- matrix(0, Dd, D)
sens_mis_se <- spec_mis_se <- matrix(0, Dd, D)
mis_pr <- mis_pr_se <- matrix(0, Dd, D)
n_vals <- 0:thres
p_vals <- (thres + 1):total_count
for (i in 1:D) {
    pars <- design[i, ]
    gamma0 <- pars[1]
    beta0 <- pars[2]
    beta_R <- pars[3]
    X <- matrix(rnorm(N * P), N, P) / sqrt(P)
    mp <- 1 / (1 + exp(-(gamma0 + X %*% gamma)))
    y_mean_n <- 1 / (1 + exp(-(beta0 + X %*% beta)))
    y_mean_p <- 1 / (1 + exp(-(beta0 + beta_R + X %*% beta)))
    yprob_n <- pbinom(thres, total_count, y_mean_n)
    yprob_p <- 1 - pbinom(thres, total_count, y_mean_p)
    pr[i] <- mean(mp)
    pr_se[i] <- sd(mp) / sqrt(N)
    sens[i] <- mean(yprob_p * mp) / mean(yprob_p * mp + (1 - yprob_n) * (1 - mp))
    sens_se[i] <- ratio_se(yprob_p * mp, yprob_p * mp + (1 - yprob_n) * (1 - mp))
    spec[i] <- mean(yprob_n * (1 - mp)) / mean(yprob_n * (1 - mp) + (1 - yprob_p) * mp)
    spec_se[i] <- ratio_se(yprob_n * (1 - mp), yprob_n * (1 - mp) + (1 - yprob_p) * mp)
    for (j in 1:Dd) {
        pars_mis <- design_mis[j, ]
        alpha0 <- pars_mis[1]
        alpha_Y <- pars_mis[2]
        alpha_R <- pars_mis[3]
        # Need to marginalize cond'l on the fact Y <= thres / Y > thres
        o_nn <- rowSums(sapply(n_vals, function(y) dbinom(y, total_count, y_mean_n) / (1 + exp(-(alpha0 + y * alpha_Y + X %*% alpha))))) # nolint
        on_nn <- rowSums(sapply(n_vals, function(y) dbinom(y, total_count, y_mean_n) / (1 + exp(alpha0 + y * alpha_Y + X %*% alpha)))) # nolint
        o_np <- rowSums(sapply(n_vals, function(y) dbinom(y, total_count, y_mean_p) / (1 + exp(-(alpha0 + y * alpha_Y + alpha_R + X %*% alpha))))) # nolint
        on_np <- rowSums(sapply(n_vals, function(y) dbinom(y, total_count, y_mean_p) / (1 + exp(alpha0 + y * alpha_Y + alpha_R + X %*% alpha)))) # nolint
        o_pn <- rowSums(sapply(p_vals, function(y) dbinom(y, total_count, y_mean_n) / (1 + exp(-(alpha0 + y * alpha_Y + X %*% alpha))))) # nolint
        on_pn <- rowSums(sapply(p_vals, function(y) dbinom(y, total_count, y_mean_n) / (1 + exp(alpha0 + y * alpha_Y + X %*% alpha)))) # nolint
        o_pp <- rowSums(sapply(p_vals, function(y) dbinom(y, total_count, y_mean_p) / (1 + exp(-(alpha0 + y * alpha_Y + alpha_R + X %*% alpha))))) # nolint
        on_pp <- rowSums(sapply(p_vals, function(y) dbinom(y, total_count, y_mean_p) / (1 + exp(alpha0 + y * alpha_Y + alpha_R + X %*% alpha)))) # nolint
        sens_obs[j, i] <- mean(o_pp * mp) / mean(o_pp * mp + o_pn * (1 - mp)) # nolint
        sens_obs_se[j, i] <- ratio_se(o_pp * mp, o_pp * mp + o_pn * (1 - mp)) # nolint
        sens_mis[j, i] <- mean(on_pp * mp) / mean(on_pp * mp + on_pn * (1 - mp)) # nolint
        sens_mis_se[j, i] <- ratio_se(on_pp * mp, on_pp * mp + on_pn * (1 - mp)) # nolint
        spec_obs[j, i] <- mean(o_nn * (1 - mp)) / mean(o_nn * (1 - mp) + o_np * mp) # nolint
        spec_obs_se[j, i] <- ratio_se(o_nn * (1 - mp), o_nn * (1 - mp) + o_nn * mp) # nolint
        spec_mis[j, i] <- mean(on_nn * (1 - mp)) / mean(on_nn * (1 - mp) + on_np * mp) # nolint
        spec_mis_se[j, i] <- ratio_se(on_nn * (1 - mp), on_nn * (1 - mp) + on_np * mp) # nolint
        mis_pr[j, i] <- 1 - mean(o_pp * mp + o_pn * (1 - mp) + o_nn * (1 - mp) + o_np * mp) # nolint
        mis_pr_se[j, i] <- sd(o_pp * mp + o_pn * (1 - mp) + o_nn * (1 - mp) + o_np * mp) # nolint
    }
}
rbind(sens, sens_se, spec, spec_se, pr, pr_se)
# rbind(mis_pr, mis_pr_se)
# rbind(
#     sens_obs[1, ], sens_obs_se[1, ],
#     sens_mis[1, ], sens_mis_se[1, ],
#     spec_obs[1, ], spec_obs_se[1, ],
#     spec_mis[1, ], spec_mis_se[1, ]
# )
# rbind(
#     sens_obs[2, ], sens_obs_se[2, ],
#     sens_mis[2, ], sens_mis_se[2, ],
#     spec_obs[2, ], spec_obs_se[2, ],
#     spec_mis[2, ], spec_mis_se[2, ]
# )
