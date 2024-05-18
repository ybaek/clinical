log_sum_exp <- function(v) {
    vm <- max(v)
    log(vm) + log(sum(exp(v - vm)))
}

inv_logit <- function(v) {
    1 / (1 + exp(-v))
}

inv_logit1m <- function(v) {
    exp(-v) / (1 + exp(-v))
}

dbinom_logit <- function(x, size, par) {
    dbinom(x, size, inv_logit(par))
}

pbinom_logit <- function(q, size, par) {
    pbinom(q, size, inv_logit(par))
}

qbinom_logit <- function(p, size, par) {
    qbinom(p, size, inv_logit(par))
}

rbinom_logit <- function(n, size, par) {
    rbinom(n, size, inv_logit(par))
}

generate_R <- function(X, gamma, gamma0) {
    logit_mean <- gamma0 + X %*% gamma
    n <- length(logit_mean)
    R <- rbinom_logit(n, 1, logit_mean)
    R
}

generate_Y <- function(X, R, total_count, gamma, gamma0, beta, beta0, beta_R) {
    logit_mean <- beta0 + X %*% beta + R * beta_R
    n <- length(logit_mean)
    Y <- rbinom_logit(n, total_count, logit_mean)
    Y
}

generate_O <- function(X, Y, R, alpha, alpha0, alpha_Y, alpha_R, alpha_YR) {
    logit_mean <- alpha0 + X %*% alpha + Y * alpha_Y + R * alpha_R +
        Y * R * alpha_YR
    n <- length(logit_mean)
    O <- rbinom_logit(n, 1, logit_mean)
    O
}

ratio_se <- function(num, denom) {
    # Delta method approximation for ratio of Monte Carlo estimators
    sqrt(
        var(num) / mean(denom)^2 -
            2 * mean(num) / mean(denom)^3 * cov(num, denom) +
            var(denom) * mean(num)^2 / mean(denom)^4
    ) / sqrt(length(num))
}
