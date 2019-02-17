function THETAmatrix = pmcmc_sdmem_tumor(problem,group,data,startstate,starttime,bigtheta,parmask,parbase,R_mcmc,step_rw,covariance_refresh,numparticles,numAPFsim)

% a particle marginal method (denoted PMCMC in the paper) using an auxiliary particle filter sequential Monte Carlo method to obtain an
% unbiased approximation of the likelihood function.
%
% Input:
%      - problem: a string, see tumor_run.m
%      - group: useless, as we should only pass data cotaining exclsively the corresponding group to fit
%      - data: data for the given group, arranged as a n x 3 matrix, where column 1 has sampling times for all subject in the chosen group,
%              column 2 has corresponding log-volumes, and column 3 has a
%              subjects ID.
%      - startstate: vector of initial volums v0, for each subject
%      - starttime: starting time to initialize the simulations, e.g time=0
%      - bigtheta: vector with all model parameters (both free and fixed ones)
%      - parmask: a vector of zeroes and ones. Zeroes denote fixed parameters that need no estimation (these have values set in
%                 parbase), ones denote free to vary parameters that are object of inference.
%      - parbase: safe to set as equal to bigtheta, may contain fixed constants that we need not to estimate (hence these have
%                 corresponding entries = 0 in parmask)
%      - R_mcmc: the number of MCMC iterations.
%      - step_rw: initial values for the standard deviations of the adaptive Gaussian random walk (Haario et al., Bernoulli 2001)
%      - covariance_refresh: how often we update the covariance matrix in adaptive Metropolis  random walk.
%      - numparticles: number of particles for the auxiliary particle filter (denoted L in the paper)
%      - numAPFsim: number of particles required to compute the means in the first stage of the auxiliary particle filter (denoted L2 in the paper)
%Output:
%      - THETAmatrix: a R_mcmc x sum(parmask) matrix of posterior draws, one column for each parameter to be inferred


fprintf('\nSimulation has started...')
% extract parameters to be estimated (i.e. those having parmask==1, from the full vector bigtheta)
theta_old = param_mask(bigtheta,parmask);
MCMC = zeros(R_mcmc,length(theta_old)); 
MCMC(1,:) = theta_old;

%:::::::::::::::::::::::::::::: INITIALIZATION  ::::::::::::::::

bigtheta_old = bigtheta;

loglik_old = tumor_smclikelihood(startstate,starttime,bigtheta_old,data,numparticles,numAPFsim);

%:::::::::::::::::::::::::::::: PRIORS :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   % evaluate priors at old parameters
   prior_old =  feval([problem, '_prior'],theta_old);
   % old priors product
   prod_priors_old = prior_old;  % this assignment is redundant...kept for compatibility with older versions
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if ~isfinite(loglik_old) || prod_priors_old==0 
    error('The initial proposal is not admissible.')
end


% initial (diagonal) covariance matrix for the Gaussian proposal
cov_current = diag(step_rw.^2);
% propose a value for parameters using Gaussian random walk
theta = mvnrnd(theta_old,cov_current);

% reconstruct updated full vector of parameters (i.e. rejoins
% fixed/constant parameters and those to be estimated)
bigtheta = param_unmask(theta,parmask,parbase);
loglik = tumor_smclikelihood(startstate,starttime,bigtheta,data,numparticles,numAPFsim);

%:::::::::::::::::::::::::::::: PRIORS :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   % evaluate priors at proposed parameters
   prior = feval([problem, '_prior'],theta);
   % proposal priors product
   prod_priors = prior ;  % this assignment is redundant...kept for compatibility with older versions
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 if log(rand) < loglik-loglik_old +log(prod_priors)-log(prod_priors_old)
      % here we accept our proposal theta
      MCMC(2,:) = theta;
      loglik_old = loglik;
      theta_old = theta;
      prod_priors_old = prod_priors;
 else
     % reject proposal
      MCMC(2,:) = theta_old;
 end



accept_proposal=0;  % start the counter for the number of accepted proposals
num_proposal=0;     % start the counter for the total number of proposed values

length_CoVupdate = covariance_refresh;  % the frequency of covariance update in adaptive Metropolis


for mcmc_iter = 3:R_mcmc
   
    %::::::::: ADAPTATION OF THE COVARIANCE MATRIX FOR THE PARAMETERS PROPOSAL :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    %::::::::: here we follow the adaptive Metropolis method as in:
    %::::::::  Haario et al. (2001) "An adaptive Metropolis algorithm", Bernoulli Volume 7, 223-242.



       if mcmc_iter == length_CoVupdate 
             lastCovupdate = 0;
             cov_last = cov_current;
             cov_old = cov_last;
             sum_old = theta_old;
       end
       if (mcmc_iter < length_CoVupdate)
          cov_current = diag(step_rw.^2); 
          theta = mvnrnd(theta_old,cov_current);
       else
           if (mcmc_iter == lastCovupdate+length_CoVupdate) 
               % we do not need to recompute the covariance on the whole
               % past history in a brutal way: we can use a recursive
               % formula. See the reference in the file cov_update.m
               covupdate = cov_update(MCMC(lastCovupdate+1:mcmc_iter-1,:),sum_old,length_CoVupdate-1,cov_old);
               sum_old = sum(MCMC(lastCovupdate+1:mcmc_iter-1,:));
               cov_old = cov(MCMC(lastCovupdate+1:mcmc_iter-1,:));
               % compute equation (1) in Haario et al.
               cov_current = (2.38^2)/length(theta)*covupdate +  (2.38^2)/length(theta) * 1e-12 * eye(length(theta)); 
               theta = mvnrnd(theta_old,cov_current);
               cov_last = cov_current;
               lastCovupdate = mcmc_iter;
               fprintf('\nPMCMC iteration -- adapting covariance...')
               fprintf('\nPMCMC iteration #%d -- acceptance ratio %4.3f percent',mcmc_iter,accept_proposal/num_proposal*100)
               accept_proposal=0;
               num_proposal=0;
               THETAmatrix_temp = MCMC(1:mcmc_iter-1,:);
               save('THETAmatrix_temp','THETAmatrix_temp') 
           else
              % Here there is no "adaptation" for the covariance matrix,
              % hence we use the same one obtained at last update
                theta = mvnrnd(theta_old,cov_last);
           end
       end
    %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    %::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::   

    num_proposal = num_proposal+1;
    
    bigtheta = param_unmask(theta,parmask,parbase);
    loglik = tumor_smclikelihood(startstate,starttime,bigtheta,data,numparticles,numAPFsim);

    % evaluate priors at proposed parameters
    prior = feval([problem, '_prior'],theta);
    prod_priors = prior;
    
      %   logkernel-logkernel_old +log(prod_priors)-log(prod_priors_old)
    if log(rand) < loglik-loglik_old +log(prod_priors)-log(prod_priors_old)
             accept_proposal=accept_proposal+1;
             MCMC(mcmc_iter,:) = theta;
             loglik_old = loglik;
             theta_old = theta;
             prod_priors_old = prod_priors;
    else
             MCMC(mcmc_iter,:) = theta_old;
    end
   
end

THETAmatrix = MCMC;
   
save('THETAmatrix','THETAmatrix')    
    
end



