data { 
  int P; // covariate dimension (incl. R at the end)
  int T; // total binomial count
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  int<lower=1,upper=N_obs+N_mis> ii_obs[N_obs];
  int<lower=1,upper=N_obs+N_mis> ii_mis[N_mis];
  int<lower=0,upper=1> r[N_obs+N_mis];
  matrix[N_obs+N_mis,P] x;
  int<lower=1,upper=T> y_obs[N_obs];
}
transformed data {
  int<lower=0> N = N_obs+N_mis;
  matrix[N,P+1] xr = append_col(x, to_vector(r));
  int<lower=0,upper=1> o[N];
  o[ii_obs] = rep_array(1, N_obs);
  o[ii_mis] = rep_array(0, N_mis);
}
parameters {
  // Model r ~ x
  real gamma0;
  vector[P] gamma;
  // Model y ~ r + x
  real beta0;
  real beta_r;
  vector[P] beta;
}
transformed parameters {
  vector[P+1] beta_c = append_row(beta, beta_r);
}
model { 
  // constants
  array[N] int ones_r = ones_int_array(N);
  array[N_obs] int Ts_y1 = rep_array(T, N_obs);
  array[N_obs] int ones_o = ones_int_array(N_obs);
  // Prior
  gamma0 ~ std_normal();
  gamma ~ std_normal();
  beta0 ~ std_normal();
  beta_r ~ std_normal();
  beta ~ std_normal();
  // Likelihood
  target += binomial_logit_lpmf(r | ones_r, gamma0 + x*gamma) +
    binomial_logit_lpmf(y_obs | Ts_y1, beta0 + xr[ii_obs,:]*beta_c);
}
