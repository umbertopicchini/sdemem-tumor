
data {
  int<lower=0> N;
  int<lower=0> M;
  real x[M];
  int<lower=0> t[N+1];
  real y[M];
  real logv0[N];  
}
parameters {
  real<lower=0, upper=1> mu_alpha;
  real mu_logbeta;    
  real mu_logdelta;            

  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_delta;
  real<lower=0> sigma_epsilon;

  real<lower=0, upper=1> alpha[N];  
  real <lower=0> beta[N];
  real <lower=0> delta[N];

}
transformed parameters {
    real mu_beta;
    real mu_delta;

    mu_beta = exp(mu_logbeta);
    mu_delta = exp(mu_logdelta);
}
model {
  mu_alpha ~ normal(0.6, 0.2)T[0, 1]; 
  mu_logbeta ~ normal(0.7, 0.6);  
  mu_logdelta ~ normal(0.7, 0.6); 
  sigma_alpha ~ inv_gamma(5, 1.5);  
  sigma_beta ~ inv_gamma(4, 2);   
  sigma_delta ~ inv_gamma(4, 2);   
  sigma_epsilon ~ inv_gamma(2, 1);  
  for (n in 1:N)
     alpha[n] ~ normal(mu_alpha, sigma_alpha)T[0, 1]; 
  beta ~ normal(mu_beta, sigma_beta);  
  delta ~ normal(mu_delta, sigma_delta);  
  
  for (n in 1:N)
    for (j in t[n]+1:t[n+1]) 
        y[j] ~ normal(logv0[n]+log((1-alpha[n])*exp(beta[n]*x[j])+alpha[n]*exp(-1*delta[n]*x[j])), sigma_epsilon);

}
