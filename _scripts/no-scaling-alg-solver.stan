//'
//' Ill-fated (at least for now) using stan's algebraic solver for this problem.
//'
functions {
    vector beta_starts(real shape,
                       real offset,
                       real total_N0,
                       int n_stages) {

        vector times[n_stages+1];
        times[1] = offset;
        real t_step = 1.0 / (n_stages * 1.0);  // multiplying by 1.0 to convert to real
        for (i in 2:(n_stages+1)) times[i] = times[i-1] + t_step;

        vector aphids0[n_stages];
        real sum_aphids0 = 0;

        real corr_time, pbeta_val, pbeta_val0;
        for (i in 1:(n_stages+1)) {
            corr_time = times[i];
            if (times[i] > 1) corr_time -= 1.0;
            pbeta_val = beta_cdf(corr_time, shape, shape);
            if (times[i] > 1) pbeta_val += 1.0;
            if (i > 1) {
                aphids0[i-1] = pbeta_val - pbeta_val0;
                sum_aphids0 += aphids0[i-1];
            }
            pbeta_val0 = pbeta_val;
        }

        aphids0 /= sum_aphids0;
        aphids0 *= total_N0;

        return aphids0;
    }

    vector sim_pcg(vector abunds0,
                   array[] int time,
                   matrix L,
                   real K) {
        int n_obs = num_elements(time);
        vector abunds = abunds0;
        real N_t = sum(abunds0);    // N at time t
        real N_t1 = N_t;            // N at time t+1
        vector pcg[n_obs];
        int max_t = max(time);
        real S;
        int i = 1;
        for (t in 1:max_t) {
            S = 1 / (1 + N_t / K);
            abunds = S * (L * abunds);
            N_t1 = sum(abunds);
            if (t == time[i]) {
                pcg[i] = N_t1 / N_t;
                i++;
            }
            N_t = N_t1;
        }
        return pcg;
    }

    vector system(vector y,              // unknowns
                  vector theta,          // parameters
                  data array[] real x_r, // data (real)
                  array[] int x_i) {     // data (integer)
        real shape = theta[1];
        real offset = theta[2];
        vector[2] z;
        z[1] = y[1] - theta[1];
        z[2] = y[1] * y[2] - theta[2];
        return z;
    }
}
data {
    // indices
    int<lower=1> n_obs;             // # observations
    int<lower=1> n_stages;          // # stages
    real<lower=1> total_N0;         // starting density across all stages
    real<lower=1> K;                // density dependence
    // data
    matrix<lower=0> L[n_stages, n_stages];  // Leslie matrix
    vector pcg[n_obs];                      // per-capita growth rates
    int<lower=1> time[n_obs];               // times
}
parameters {
    // real alpha[n_coef];                         // fixed effects and intercepts
    // real z[sum(lev_per_g)];                     // standardized variates for group levels
    // real<lower=0, upper=1> phi[max(p_groups)];  // autoregressive parameter for each
    // real<lower=0> sig_beta[sum(g_per_ff)];      // group standard deviations
    // real<lower=0> sig_res;                      // residual standard deviation
    real<lower=1> shape;
    real<lower=0, upper=1> offset;
    real<lower=0> sigma;
}
transformed parameters {
    vector<lower=0> abunds0[n_stages];      // abundance by stage to start
    abunds0 = beta_starts(shape, offset, total_N0, n_stages);
    vector y_pred[n_obs];     // predicted values
    y_pred = sim_pcg(abunds0, time, L, K);
}
model {
    // // priors:
    // alpha ~ normal(0, 1);
    // z ~ normal(0, 1);
    // for (i in 1:sum(g_per_ff)){
    //     sig_beta[i] ~ gamma(1.5, 3);
    // }
    // for (i in 1:max(p_groups)){
    //     // phi[i] ~ normal(0, 0.5) T[0, p_bound];
    //     phi[i] ~ beta(2, 2);
    // }
    // sig_res ~ gamma(1.5, 3);
    // // observations:
    // y ~ normal(y_pred, sig_res);
    shape ~ gamma(1.5, 3) T[1.0,];
    offset ~ beta(2, 2);
    sigma ~ gamma(1.5, 3);
    // observations:
    y ~ normal(y_pred, sigma);
}
generated quantities {
  real log_lik = 0;
  for(i in 1:n_obs){
    log_lik += normal_lpdf(y[i] | y_pred[i], sigma);
  }
}
