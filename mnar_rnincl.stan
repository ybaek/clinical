data { 
  int P; // covariate dimension (incl. R at the end)
  int T; // total binomial count
  int thres; // thresholding y
  int<lower=0> N_obs;
  int<lower=0> N_mis;
  int<lower=1,upper=N_obs+N_mis> ii_obs[N_obs];
  int<lower=1,upper=N_obs+N_mis> ii_mis[N_mis];
  matrix[N_obs+N_mis,P] x;
  int<lower=1,upper=T> y_obs[N_obs];
}
transformed data {
  int<lower=0> N = N_obs+N_mis;
  int<lower=0,upper=1> z_obs[N_obs];
  for (n in 1:N_obs) {
    z_obs[n] = y_obs[n] > thres;
  }
  matrix[N_obs,P+1] xz = append_col(x[ii_obs,:], to_vector(z_obs));
  int<lower=0,upper=1> o[N];
  o[ii_obs] = rep_array(1, N_obs);
  o[ii_mis] = rep_array(0, N_mis);
}
parameters {
  // Model y ~ x
  real beta0;
  vector[P] beta;
  // Model obs/mis ~ z + x
  real alpha0;
  real alpha_z;
  vector[P] alpha;
}
transformed parameters {
  vector[P+1] alpha_c = append_row(alpha, alpha_z);
}
model { 
  // constants
  array[N] int ones_r = ones_int_array(N);
  array[N_obs] int Ts_y1 = rep_array(T, N_obs);
  array[N_mis] int Ts_y2 = rep_array(T, N_mis);
  array[N_obs] int ones_o = ones_int_array(N_obs);
  array[N_mis] int ones_m = ones_int_array(N_mis);
  // Prior
  beta0 ~ std_normal();
  beta ~ std_normal();
  alpha0 ~ std_normal();
  alpha_z ~ std_normal();
  alpha ~ std_normal();
  // Likelihood
  target += binomial_logit_lpmf(y_obs | Ts_y1, beta0 + x[ii_obs,:]*beta) +
    binomial_logit_lpmf(o[ii_obs] | ones_o, alpha0 + xz*alpha_c);
  vector[T] ll_ymis;
  for (k in 1:T) {
    array[N_mis] int ks = rep_array(k, N_mis);
    int zk = k > thres;
    real intercept_k = alpha0 + zk*alpha_z;
    ll_ymis[k] = binomial_logit_lpmf(o[ii_mis] | ones_m, intercept_k + x[ii_mis,:]*alpha) +
      binomial_logit_lpmf(ks | Ts_y2, beta0 + x[ii_mis,:]*beta);
  }
  target += log_sum_exp(ll_ymis);
}
