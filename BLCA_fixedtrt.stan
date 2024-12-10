data {
    int<lower = 1> N;                   // number of individuals
    int<lower = 1> K;                   // number of covariates
    int<lower = 1> C;                   // number of latent classes
    matrix[N, K] x;                     // covariates
    vector[N] y;                        // outcomes
    int<lower = 0, upper = 1> t[N];     // treatment indicators
}

transformed data { // standardize the data
    vector[K] mu_x;
    vector[K] sd_x;
    matrix[K, N] x_std;
    real mu_y = mean(y);
    real sd_y = sd(y);
    vector[N] y_std = (y - mu_y) / sd_y;

    x_std[1, ] = to_row_vector(x[, 1]);
    for (k in 2:K) {
        mu_x[k] = mean(x[, k]);
        sd_x[k] = sd(x[, k]);
        x_std[k, ] = (to_row_vector(x[, k]) - mu_x[k]) / sd_x[k];
    }
}

parameters {
    matrix[K, C - 1] beta_c_std;          // mixing probability model
    matrix[K, C] beta_t_std;              // selection model - not needed at this time due to random treatment assignment
    matrix[K, C] beta_y_std;              // outcome model
    vector[C] tau_std;                    // average treatment effects
    vector<lower = 0>[C] sigma_y_std;     // standard deviation of error term
}

transformed parameters {
    vector[N] log_lik;
    vector[C] pr[N];                      // probabilities of membership
    vector[C] delta[N];                   // doubly robust estimator for each unit
    vector[C] dr;                         // doubly robust estimator
    vector[C] dr_trim;                    // doubly robust estimator (trimmed)

    for (n in 1:N) {
        vector[C] q;
        q[1] = 0;
        for (c in 2:C) {
            q[c] = dot_product(x_std[, n], beta_c_std[, c - 1]);
        }
        q = log(softmax(q));

        for (c in 1:C) {
            real p = 0.5; // Fixed treatment assignment probability
            real y0_std = dot_product(x_std[, n], beta_y_std[, c]);
            real y1_std = y0_std + tau_std[c];
            real delta_std = (t[n] / p * y_std[n] + (p - t[n]) / p * y1_std)
                            - ((1 - t[n]) / (1 - p) * y_std[n] + (t[n] - p) / (1 - p) * y0_std);
            // Remove ps[n, c] as it is not needed anymore
            q[c] += bernoulli_lpmf(t[n] | p) + normal_lpdf(y_std[n] | y0_std + t[n] * tau_std[c], sigma_y_std[c]);
            delta[n, c] = sd_y * delta_std;
        }
        log_lik[n] = log_sum_exp(q);
        pr[n] = exp(q - log_lik[n]);
    }

    // Doubly robust estimator is not strictly needed here
    // If you wish to keep it, it should be adjusted as follows:
    for (c in 1:C) {
    real a = 0;
    real b = 0;
    real a_trim = 0;
    real b_trim = 0;
    for (n in 1:N) {
        a += delta[n, c] * pr[n, c];
        b += pr[n, c];
        if (0.1 < delta[n, c] && delta[n, c] < 0.9) {
            a_trim += delta[n, c] * pr[n, c];
            b_trim += pr[n, c];
        }
    }
    dr[c] = a / b;
    dr_trim[c] = b_trim > 0 ? a_trim / b_trim : 0;  // Safeguard against division by zero
}
}


// Define priors
model {
    to_vector(beta_c_std) ~ normal(0, 1); 
    to_vector(beta_t_std) ~ normal(0, 1);  // Not used but kept for consistency
    to_vector(beta_y_std) ~ normal(0, 1);
    tau_std ~ normal(0, 1);
    sigma_y_std ~ normal(0, 1) T[0, ]; // cauchy(0, 1);
    target += sum(log_lik);
}

generated quantities {
    matrix[K, C - 1] beta_c;
    matrix[K, C] beta_t;
    matrix[K, C] beta_y;
    vector[C] tau = tau_std * sd_y;
    vector[C] sigma_y = sigma_y_std * sd_y;

    beta_c[1, ] = beta_c_std[1, ];
    beta_t[1, ] = beta_t_std[1, ];
    beta_y[1, ] = sd_y * beta_y_std[1, ] + mu_y;

    for (k in 2:K) {
        beta_c[k, ] = beta_c_std[k, ] / sd_x[k];
        beta_c[1, ] -= beta_c[k, ] * mu_x[k];
        beta_t[k, ] = beta_t_std[k, ] / sd_x[k];
        beta_t[1, ] -= beta_t[k, ] * mu_x[k];
        beta_y[k, ] = sd_y / sd_x[k] * beta_y_std[k, ];
        beta_y[1, ] -= beta_y[k, ] * mu_x[k];
    }
}


