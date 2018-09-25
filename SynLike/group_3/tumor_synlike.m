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
logdelta = bigtheta(2);
logalpha = bigtheta(3);
loggamma = bigtheta(4);
logtau = bigtheta(5);
logsigmabeta = bigtheta(6);
logsigmadelta = bigtheta(7);
logsigmaalpha = bigtheta(8);
logsigmaerror = bigtheta(9);

% exponentiate, for log-parameters
beta = exp(logbeta);
delta = exp(logdelta);
alpha = exp(logalpha);
gamma = exp(loggamma);
tau = exp(logtau);
sigmabeta = exp(logsigmabeta);
sigmadelta = exp(logsigmadelta);
sigmaalpha = exp(logsigmaalpha);
sigmaerror = exp(logsigmaerror);



obssummaries = tumor_summaries(alllogdata(:,2),alllogdata(:,1),alllogdata(:,3));  % observed summaries
dsum = length(obssummaries);

lognoisyobs = [];
allvolume = [];
subjectsid = unique(alllogdata(:,3));


count_subj=1;
for subject = subjectsid'   % NEAT! if I take the transpose of subjectsid (so now it's a row vector) I can iterate through it
    sublogdata = alllogdata(alllogdata(:,3)==subject,:); %extract data pertaining the current subject
    n = size(sublogdata,1);  % number of measurements for the given subject
   % v0 = sublogdata(1,2); % the assumed starting volume (mm^3). Will be converted to log-scale below 
    v0 = startstate(count_subj);
    count_subj = count_subj+1;
    times = sublogdata(:,1);  % observational times for the given subject
    timesrep = repmat(times,1,numsim); 
    stepsizes = diff(timesrep,1,1);
    alpharand = mytruncgaussrandraw([0,1,alpha,sigmaalpha],[1 numsim]);  %Umberto's custom made file
    betarand = beta + sigmabeta*randn(1,numsim);
    deltarand = delta + sigmadelta*randn(1,numsim);
    v0killed = alpharand.*v0;
    v0surv = (1-alpharand).*v0;
    survived = zeros(n,numsim); % intialize matrix
    killed = survived; % intialize matrix (to zeros)
    volume = survived; % intialize matrix (to zeros)
    randn_std1 = randn(n,numsim);
    randn_std2 = randn(n,numsim);
    dW1 = zeros(n,numsim);             % preallocate arrays ...
    dW2 = dW1;            
    dW1(1,:) = sqrt(timesrep(1,:)-starttime).*randn_std1(1,:);      
    dW2(1,:) = sqrt(timesrep(1,:)-starttime).*randn_std2(1,:);   
    dW1(2:n,:) = sqrt(stepsizes).*randn_std1(2:n,:);      
    dW2(2:n,:) = sqrt(stepsizes).*randn_std2(2:n,:); 

    for ii=1:n
          if ii==1
              % the explicit solution of the state equations
              survived(1,:) = v0surv.*exp(betarand.*(timesrep(1,:)-starttime)+gamma*dW1(1,:));
              killed(1,:) = v0killed.*exp(-deltarand.*(timesrep(1,:)-starttime)+tau*dW2(1,:));
              volume(1,:) = survived(1,:) + killed(1,:);
          else
              % propagate forward from the particles having sampled indeces 
              survived(ii,:) = survived(ii-1,:).*exp(betarand.*stepsizes(ii-1,:)+gamma*dW1(ii,:));
              killed(ii,:) = killed(ii-1,:).*exp(-deltarand.*stepsizes(ii-1,:)+tau*dW2(ii,:));
              volume(ii,:) = survived(ii,:) + killed(ii,:);
          end
    end

   allvolume = [allvolume;volume];  % stack all volumes 

end


   % transform to log scale for consistency with data then add measurement error
   lognoisyobs = log(allvolume) + sigmaerror*randn(size(allvolume,1),numsim);

   simsummaries = tumor_summaries(lognoisyobs,alllogdata(:,1),alllogdata(:,3));  % simulated summaries
   mean_simsummaries = mean(simsummaries,2);
   %cov_summaries = (summaries_individual-repmat(mean_summaries,1,numsim))*(summaries_individual-repmat(mean_summaries,1,numsim))' / (numsim-1);
   cov_simsummaries = cov(simsummaries');

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


