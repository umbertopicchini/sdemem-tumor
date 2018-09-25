
data {
  int<lower=0> N;
  int<lower=0> M;
  real x[M];
  int<lower=0> t[N+1];
  real y[M];
  real logv0[N];  // THESE ARE NOW GIVEN AS DATA FROM THE RUN FILE
}
parameters {
  real<lower=0, upper=1> mu_alpha;
  real mu_logbeta;    // THE MEAN OF THE \LOGBETA PARAMETER      
  real mu_logdelta;   // THE MEAN OF THE \LOGDELTA PARAMETER         
//  real mu_v0;          

//  real<lower=0> sigma_v0;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_delta;
  real<lower=0> sigma_epsilon;

  real<lower=0, upper=1> alpha[N];  
  real <lower=0> beta[N];
  real <lower=0> delta[N];
 // real logv0[N];

}
transformed parameters {
//  real<lower=0> sigma_v0;       // WE DO NOT TRANSFORM THOSE ANY LONGER
//  real<lower=0> sigma_alpha;
//  real<lower=0> sigma_beta;
//  real<lower=0> sigma_delta;
//  real<lower=0> sigma_epsilon;
//  real b0[N];
//  real d0[N];
//  real a2[N,2];
//  real a1[N,2];
//  real a0[N];
    real mu_beta;
    real mu_delta;

    mu_beta = exp(mu_logbeta);
    mu_delta = exp(mu_logdelta);
 // sigma_v0 = sqrt(sigmasq_v0);
 // sigma_alpha = sqrt(sigmasq_alpha);
//  sigma_beta = sqrt(sigmasq_beta);
//  sigma_delta = sqrt(sigmasq_delta);
//  sigma_epsilon = sqrt(sigmasq_epsilon);
//  for (n in 1:N)
//  a0[n] = exp(logalpha[n]);
//  for (n in 1:N)
//  b0[n] = exp(logbeta[n]);
//  for (n in 1:N)
//  d0[n] = exp(logdelta[n]);
 // for (n in 1:N)
//  a2[n,1] = exp(logalpha[n]);  
//  for (n in 1:N)
//  a2[n,2] = 1.0;
//  for (n in 1:N)
//  a1[n,1] = min(a2[n,]);
//  for (n in 1:N)
//  a1[n,2] = 0.0;
//  for (n in 1:N)
//  a0[n] = max(a1[n,]);
}
model {
//  mu_v0 ~ normal(5, 1);
  mu_alpha ~ normal(0.6, 0.2)T[0, 1];  // THE SECOND ARGUMENT SHOULD BE A SD NOT A VARIANCE. AND SHOULD BE TRUNCATED TO [0,1]
  mu_logbeta ~ normal(0.7, 0.6);  // THE SECOND ARGUMENT SHOULD BE A SD NOT A VARIANCE
  mu_logdelta ~ normal(0.7, 0.6); // THE SECOND ARGUMENT SHOULD BE A SD NOT A VARIANCE
//  sigma_v0 ~ inv_gamma(4, 2);
  sigma_alpha ~ inv_gamma(5, 1.5);  // IN THE PAPER WE HAVE INVGAMMA PRIORS FOR STD NOT VARIANCES
  sigma_beta ~ inv_gamma(4, 2);    // IN THE PAPER WE HAVE INVGAMMA PRIORS FOR STD NOT VARIANCES
  sigma_delta ~ inv_gamma(4, 2);   // IN THE PAPER WE HAVE INVGAMMA PRIORS FOR STD NOT VARIANCES
  sigma_epsilon ~ inv_gamma(2, 1);  // IN THE PAPER WE HAVE INVGAMMA PRIORS FOR STD NOT VARIANCES
//  logv0 ~ normal(mu_v0, sigma_v0);
  for (n in 1:N)
     alpha[n] ~ normal(mu_alpha, sigma_alpha)T[0, 1]; // SHOULD BE TRUNCATED ON [0,1]
  beta ~ normal(mu_beta, sigma_beta);  // vectorized
  delta ~ normal(mu_delta, sigma_delta);  // vectorized
  
  for (n in 1:N)
    for (j in t[n]+1:t[n+1]) 
    //  y[j] ~ normal(logv0[n]+log(a0[n]*exp(b0[n]*x[j])+(1-a0[n])*exp(-1*d0[n]*x[j])), sigma_epsilon);
    // SHOULD BECOME
        y[j] ~ normal(logv0[n]+log((1-alpha[n])*exp(beta[n]*x[j])+alpha[n]*exp(-1*delta[n]*x[j])), sigma_epsilon);

}
