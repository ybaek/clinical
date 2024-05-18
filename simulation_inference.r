source("prior_model.r")
# Simulation DGM
P <- 30
total_count <- 40
thres <- 14
N <- 300
set.seed(202453001)
X <- matrix(rnorm(N * P), N) / sqrt(P)
par_grid <- data.frame(
    beta0   = rep(c(1, -0.5, 0), 2),
    beta_r  = rep(c(0.2, 0.05, 0.01), 2),
    gamma0  = rep(c(1.5, 1.5, -1.4), 2),
    alpha0  = rep(0, 6),
    alpha_z = rep(c(0.1, -1.2), each = 3)
)
par_grid$setting <- seq(dim(par_grid)[1])
beta <- rnorm(P) / sqrt(P)
gamma <- sqrt(.8) * beta + rnorm(P, 0, sqrt(.2 / P))
alpha <- rnorm(P) / sqrt(P)

# Fitting the model in Stan
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(rstan)
options(mc.cores = parallel::detectCores() - 1)
rstan_options(auto_write = TRUE)

parset1 <- "beta0"
parset2 <- c(parset1, c("beta_r", "gamma0"))
parset3 <- c(parset1, c("alpha0", "alpha_z"))
parset4 <- c(parset2, c("alpha0", "alpha_z"))

for (i in 1:dim(par_grid)[1]) {
    beta0 <- par_grid[i, 1]
    beta_r <- par_grid[i, 2]
    gamma0 <- par_grid[i, 3]
    alpha0 <- par_grid[i, 4]
    alpha_z <- par_grid[i, 5]
    # 1. Generate R
    R <- generate_R(X, gamma, gamma0)
    # 2. Generate Y (total)
    Y <- generate_Y(X, R, total_count, gamma, gamma0, beta, beta0, beta_r)
    Z <- 0 + (Y > thres)
    # 3. Generate O
    O <- generate_O(X, Z, R, alpha, alpha0, alpha_z, 0, 0)
    ii_obs <- which(!!O)
    ii_mis <- which(!O)
    N_obs <- length(ii_obs)
    N_mis <- length(ii_mis)
    Yobs <- Y[ii_obs]
    data <- list(
        P = P,
        T = total_count,
        thres = thres,
        N_obs = N_obs,
        N_mis = N_mis,
        ii_obs = ii_obs,
        ii_mis = ii_mis,
        r = R,
        x = X,
        y_obs = Yobs
    )
    fit_mar <- stan(
        file = "mar_rnincl.stan",
        data = data,
        chains = 4,
        warmup = 1000,
        iter = 2000,
        pars = parset1
    )
    fit_mar_r <- stan(
        file = "mar_rincl.stan",
        data = data,
        chains = 4,
        warmup = 1000,
        iter = 2000,
        pars = parset2
    )
    fit_mnar <- stan(
        file = "mnar_rnincl.stan",
        data = data,
        chains = 4,
        warmup = 1000,
        iter = 2000,
        pars = parset3,
        init =
        )
    fit_mnar_r <- stan(
        file = "mnar_rincl.stan",
        data = data,
        chains = 4,
        warmup = 1000,
        iter = 2000,
        pars = parset4
    )
    p1 <- plot(fit_mar, pars = parset1)
    p2 <- plot(fit_mar_r, pars = parset2)
    p3 <- plot(fit_mnar, pars = parset3)
    p4 <- plot(fit_mnar_r, pars = parset4)
    tbl <- par_grid %>%
        filter(setting == i) %>%
        gather()
    p1 <- p1 +
        geom_point(
            data = tbl %>% filter(key %in% parset1),
            aes(x = value, y = rev(seq(length(parset1)))),
            size = 3, shape = 17, col = "blue"
        ) +
        scale_x_continuous(limits = c(-2, 2)) +
        ggtitle("MAR (no R)")
    p2 <- p2 +
        geom_point(
            data = tbl %>% filter(key %in% parset2),
            aes(x = value, y = rev(seq(length(parset2)))),
            size = 3, shape = 17, col = "blue"
        ) +
        scale_x_continuous(limits = c(-2, 2)) +
        ggtitle("MAR (with R)")
    p3 <- p3 +
        geom_point(
            data = tbl %>% filter(key %in% parset3),
            aes(x = value, y = rev(seq(length(parset3)))),
            size = 3, shape = 17, col = "blue"
        ) +
        scale_x_continuous(limits = c(-2, 2)) +
        ggtitle("MNAR (no R)")
    p4 <- p4 +
        geom_point(
            data = tbl %>% filter(key %in% parset4),
            aes(x = value, y = rev(seq(length(parset4)))),
            size = 3, shape = 17, col = "blue"
        ) +
        scale_x_continuous(limits = c(-2, 2)) +
        ggtitle("MNAR (with R)")
    g <- grid.arrange(p1, p2, p3, p4, nrow = 2)
    filename <- paste0("figures/fitplot", i, ".jpeg")
    ggsave(filename, g, device = "jpeg")
}
