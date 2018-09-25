function loglik = tumor_smclikelihood(startstate,starttime,bigtheta,alllogdata,numparticles,numAPFsim)
% Returns the loglikelihod corresponding to the unbiased likelihood approximation produced by an auxiliary particle filter
%
%Input:
%      - startstate: vector of initial volums v0, for each subject
%      - starttime: starting time to initialize the simulations, e.g time=0
%      - bigtheta: vector with all model parameters (both free and fixed ones)
%      - alllogdata: supplied data, see tumor_run.m
%      - numparticles: number of particles for the auxiliary particle filter (denoted L in the paper)
%      - numAPFsim: number of particles required to compute the means in the first stage of the auxiliary particle filter (denoted L2 in the paper)
%Output:
%      - an approximate loglikelihood

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

%replicas= size(timerep,2);  % timerep is a matrix that contains unique times (we assumed a balanced design --> same sampling times for all subjects), repeated for "numsubjects*replicas" number of times along the second dimension


% generate random parameters

% the function below creates a vector of numsubjects*replicas draws from a
% Gaussian with mean alpha and SD sigmaalpha truncated to [0,1]
% alpharand = alpha + TruncatedGaussian(sigmaalpha,[0,1]-alpha,[1,numsubjects_star_replicas]); % see comment from 10 April '10 in http://se.mathworks.com/matlabcentral/fileexchange/23832-truncated-gaussian

% alpharand = zeros(1,numsubjects_star_replicas);
% for ii = 1:numsubjects_star_replicas
%      alpharand(ii) = rtnorm(0,1,alpha,sigmaalpha);  % a truncated gaussian
% end


subjectsid = unique(alllogdata(:,3));
loglik = 0;  % initialize loglikelihood

count_subj = 1;
for subject = subjectsid'   % NEAT! if I take the transpose of subjectsid (so now it's a row vector) I can iterate through it
    sublogdata = alllogdata(alllogdata(:,3)==subject,:); %extract data pertaining the current subject
    n = size(sublogdata,1);  % number of measurements for the given subject
    %v0 = startstate; % the assumed starting volume (mm^3). Will be converted to log-scale below 
    v0 = startstate(count_subj);
    count_subj = count_subj+1;
    times = sublogdata(:,1);  % observational times for the given subject
    obslogvolumes = sublogdata(:,2); % measured (observed) log-volumes for the given subject
    timesrep = repmat(times,1,numparticles); 
    stepsizes = diff(timesrep,1,1);
    betarand = beta + sigmabeta*randn(1,numparticles);
    volume = zeros(n,numparticles); % intialize matrix
    particles_means = zeros(1,numparticles);
    weights_normal = 1/numparticles*ones(1,numparticles);    

    for ii=1:n
          if ii==1
              for jj=1:numparticles
                 % for each particle, propagate a further cloud of
                 % particles (this is because we want to compute central
                 % values for the auxiliary particle filter)
                 particles_volume = v0.*exp(betarand(jj).*(timesrep(1,jj)-starttime)+gamma*sqrt(timesrep(1,jj)-starttime).*randn(1,numAPFsim)); 
                 particles_means(jj) = sum(log(particles_volume))/numparticles;
              end
              % Compute importance weights (log-scale for numerical stability)
              log_firsstageweights = -1/2*log(2*pi)-log(sigmaerror) - 1/(2*sigmaerror^2)*(obslogvolumes(1) - particles_means).^2;  
              firsstageweights = exp(log_firsstageweights).*weights_normal;
              firsstageweights_norm = firsstageweights/sum(firsstageweights);
              ind_firststage = stratresample(firsstageweights_norm,numparticles); % stratified resampling of particle indeces
              volume(1,:) = v0.*exp(betarand.*(timesrep(1,:)-starttime)+gamma*sqrt(timesrep(1,:)-starttime).*randn(1,numparticles)); 
          else
              for jj=1:numparticles
                 % for each particle, propagate a further cloud of
                 % particles (this is because we want to compute central
                 % values for the auxiliary particle filter)
                 particles_volume = volume(ii-1,jj).*exp(betarand(jj).*stepsizes(ii-1,jj)+gamma*sqrt(stepsizes(ii-1,jj)).*randn(1,numAPFsim)); 
                 particles_means(jj) = sum(log(particles_volume))/numparticles;
              end
              log_firsstageweights = -1/2*log(2*pi)-log(sigmaerror) - 1/(2*sigmaerror^2)*(obslogvolumes(ii) - particles_means).^2;  
              firsstageweights = exp(log_firsstageweights).*weights_normal;
              firsstageweights_norm = firsstageweights/sum(firsstageweights);
              ind_firststage = stratresample(firsstageweights_norm,numparticles); % stratified resampling of particle indeces
              volume(ii,:) = volume(ii-1,ind_firststage).*exp(betarand.*stepsizes(ii-1,:)+gamma*sqrt(stepsizes(ii-1,:)).*randn(1,numparticles)); 
          end
    % Compute importance weights (log-scale for numerical stability)
    logweights = -1/2*log(2*pi)-log(sigmaerror) - 1/(2*sigmaerror^2)*(obslogvolumes(ii) - log(volume(ii,:))).^2 - (-1/2*log(2*pi)-log(sigmaerror) - 1/(2*sigmaerror^2)*(obslogvolumes(ii) - particles_means).^2);  
    const = max(logweights); % below  we substract the maximum value for numerical stability
    weights = exp(logweights-const); % this way it is "less likely" to obtain a numerical underflow/overflow
    % Compute loglikelihood
    loglik = loglik + const + log(sum(weights)) - log(numparticles) + log(sum(firsstageweights));        
    weights_normal = weights/sum(weights); % normalized weights
    end



end

end

