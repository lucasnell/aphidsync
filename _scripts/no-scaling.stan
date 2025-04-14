// //  #include <stan/math/prim/fun/inv_inc_beta.hpp>

functions {
    vector beta_starts(real shape,
                       real offset,
                       real total_N0,
                       int n_stages) {

        vector[n_stages+1] times;
        times[1] = offset;
        real t_step = 1.0 / (n_stages * 1.0);  // multiplying by 1.0 to convert to real
        for (i in 2:(n_stages+1)) times[i] = times[i-1] + t_step;

        vector[n_stages] aphids0;
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

    real calc_med_age(real shape,
                      real offset,
                      data int n_stages) {
        vector[n_stages] abunds = beta_starts(shape, offset, 1.0, n_stages);
        int j = 1;
        while (abunds[j] <= 0.5) {
            j += 1;
            abunds[j] += abunds[j-1];
        }
        real idx = j;
        if (j > 1 && abunds[j-1] == 0.5) {
            // This is for rare occurrence where median age is exactly between
            // two age groups:
            idx += (j - 1.0);
            idx /= 2.0;
        }
        return idx;
    }

    real calc_width99(real shape) {
        real lo = inv_inc_beta(shape, shape, 0.005);
        real hi = inv_inc_beta(shape, shape, 0.995);
        real out = abs(hi - lo);
        return out;
    }

    vector sim_pcg(vector abunds0,
                   data array[] int time,
                   data matrix L,
                   data real K) {
        int n_obs = num_elements(time);
        int n_stages = num_elements(abunds0);
        vector[n_stages] abunds = abunds0;
        real N_t = sum(abunds0);    // N at time t
        real N_t1 = N_t;            // N at time t+1
        vector[n_obs] pcg;
        int max_t = max(time);
        real S;
        int i = 1;
        for (t in 1:max_t) {
            S = 1 / (1 + N_t / K);
            abunds = S * (L * abunds);
            N_t1 = sum(abunds);
            if (t == time[i]) {
                pcg[i] = N_t1 / N_t;
                i += 1;
            }
            N_t = N_t1;
        }
        return pcg;
    }
}
data {
    // indices
    int<lower=1> n_obs;             // # observations
    int<lower=1> n_stages;          // # stages
    // data
    real<lower=1> total_N0;         // starting density across all stages
    real<lower=1> K;                // density dependence
    real<lower=1> shape_mu;         // mean for prior normal distribution for shape
    real<lower=0> shape_sigma;      // stdev for prior normal distribution for shape
    matrix<lower=0>[n_stages, n_stages] L;  // Leslie matrix
    vector[n_obs] pcg;                      // per-capita growth rates
    array[n_obs] int<lower=1> time;         // times
}
parameters {
    real<lower=1> shape;
    real<lower=0, upper=1> offset;
    real<lower=0> sigma;
}
transformed parameters {
    vector<lower=0>[n_stages] abunds0;  // abundance by stage to start
    vector[n_obs] pcg_pred;             // predicted values
    real<lower=0> med_age;
    real<lower=0> width99;
    abunds0 = beta_starts(shape, offset, total_N0, n_stages);
    pcg_pred = sim_pcg(abunds0, time, L, K);
    med_age = calc_med_age(shape, offset, n_stages);
    width99 = calc_width99(shape);
}
model {
    shape ~ normal(shape_mu, shape_sigma)T[1.0,];
    offset ~ uniform(0, 1);
    sigma ~ gamma(1.5, 3);
    // observations:
    pcg ~ normal(pcg_pred, sigma);
}
generated quantities {
  real log_lik = 0;
  for(i in 1:n_obs){
    log_lik += normal_lpdf(pcg[i] | pcg_pred[i], sigma);
  }
}
