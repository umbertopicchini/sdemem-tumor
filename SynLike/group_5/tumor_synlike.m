function [loglik,simsummaries] = tumor_synlike(startstate,starttime,bigtheta,alllogdata,numsim)

% Computes the logarithm of an unbiased estimate of a synthetic likelihood.
% This is the BSL method of Price, L. F., Drovandi, C. C., Lee, A., & Nott, D. J. (2018). Bayesian synthetic likelihood. Journal of Computational and Graphical Statistics, 27(1), 1-11.
% Input:
%       - startstate: vector of initial volums v0, for each subject
%       - starttime: starting time to initialize the simulations, e.g time=0
%       - bigtheta: vector with all model parameters (both free and fixed ones)
%       - alllogdata: supplied data, see tumor_run.m
%       - numsim: number of simulated artificial datasets ued when computing the BSL approximation
% Output:
%       - loglik: the synthetic loglikelihood
%       - simsummaries: corresponding simulated summary statistics


% get the current values of parameters. Notice these might be on log-scale.
logbeta  = bigtheta(1);
loggamma = bigtheta(2);
logsigmabeta = bigtheta(3);
logsigmaerror = bigtheta(4);

% exponentiate, for log-parameters
beta = exp(logbeta);
gamma = exp(loggamma);
sigmabeta = exp(logsigmabeta);
sigmaerror = exp(logsigmaerror);


obssummaries = tumor_summaries(alllogdata(:,2),alllogdata(:,1),alllogdata(:,3));  % observed summaries
dsum = length(obssummaries);

allvolume = [];
subjectsid = unique(alllogdata(:,3));


count_subj = 1;
for subject = subjectsid'   % NEAT! if I take the transpose of subjectsid (so now it's a row vector) I can iterate through it
    sublogdata = alllogdata(alllogdata(:,3)==subject,:); %extract data pertaining the current subject
    n = size(sublogdata,1);  % number of measurements for the given subject
    v0 = startstate(count_subj);
    count_subj = count_subj+1;
    times = sublogdata(:,1);  % observational times for the given subject
    timesrep = repmat(times,1,numsim); 
    stepsizes = diff(timesrep,1,1);
    betarand = beta + sigmabeta*randn(1,numsim);
    volume = zeros(n,numsim); % intialize matrix
    randn_std1 = randn(n,numsim);
    dW1 = zeros(n,numsim);             % preallocate arrays ...       
    dW1(1,:) = sqrt(timesrep(1,:)-starttime).*randn_std1(1,:);       
    dW1(2:n,:) = sqrt(stepsizes).*randn_std1(2:n,:);      

    for ii=1:n
          if ii==1
              volume(1,:) = v0.*exp(betarand.*(timesrep(1,:)-starttime)+gamma*dW1(1,:));
          else
              volume(ii,:) = volume(ii-1,:).*exp(betarand.*stepsizes(ii-1,:)+gamma*dW1(ii,:));
          end
    end

   allvolume = [allvolume;volume];  % stack all volumes 

end


   % transform to log scale for consistency with data then add measurement error
   lognoisyobs = log(allvolume) + sigmaerror*randn(size(allvolume,1),numsim);

   simsummaries = tumor_summaries(lognoisyobs,alllogdata(:,1),alllogdata(:,3));  % simulated summaries
   mean_simsummaries = mean(simsummaries,2);
   %cov_summaries = (summaries_individual-repmat(mean_summaries,1,numsim))*(summaries_individual-repmat(mean_summaries,1,numsim))' / (numsim-1);
   cov_simsummaries = cov(simsummaries') ;
   M = (numsim-1)*cov_simsummaries;
   
   % check positive definite covariance
   [~,positive] = chol(M) ;
   if positive>0 % leave the function as the covariance is NOT positive definite
       loglik = -inf;
       return
   end
   
   % unbiased estimator for a Gaussian density, see Price et al.
   phi_argument = M - (obssummaries-mean_simsummaries)*(obssummaries-mean_simsummaries)'/(1-1/numsim) ;
   [~,positive] = chol(phi_argument);
   if positive>0  
       loglik = -inf;
       return
   end
   loglik = -(numsim-dsum-2)/2 *logdet(M,'chol') + ((numsim-dsum-3)/2) * logdet(phi_argument,'chol');
  

end


