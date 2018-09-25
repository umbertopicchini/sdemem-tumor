function postpredcheck = tumor_posteriorpredictivecheck(startstate,starttime,THETAmatrix,alllogdata,numsim)

% Posterior predictive checks when using the BSL inference method.
% Warning: use only with data from group 1 and 3.
%
% It will simulate numposteriordraws times from the posterior
% predictive distribution of the summaries used for synthetic likelihoods
% anayses, where numposteriordraws = size(THETAmatrix,1) is the number of parameters drawn from the posterior (supplied after some burnin).
% Input: 
%       - startstate: vector of initial volums v0, for each subject
%       - starttime: starting time to initialize the simulations, e.g time=0
%       - THETAmatrix: matrix of posterior draws (one for each row), obtained using some inference method, possibly discarding burnin draws 
%       - alllogdata: supplied data, see tumor_run.m
%       - numsim: number of simulated artificial datasets ued when computing the BSL approximation
% Output:
%       - postpredcheck: a numposteriordraws x (5*numsubjects+3) matrix. Each rows contain simulated synthetic summaries, corresponding to 
%                         5 individual summaries for each subject, and additional 3 summaries that are common/shared
%                         between all subjects. E.g. group 3 has 5 subjects, therefore postpredcheck has 28 columns.
%                         See also tumor_summaries.m

numcheck = size(THETAmatrix,1);


testsummaries = tumor_summaries(alllogdata(:,2),alllogdata(:,1),alllogdata(:,3));  % simulated summaries
postpredcheck = zeros(numcheck,length(testsummaries));


fprintf('\nThis might require several minutes...')
for check = 1:numcheck

allvolume = [];
subjectsid = unique(alllogdata(:,3));
% get parameters generated from the posterior distribution (these are
% passed as input)
% get the current values of parameters. Notice these might be on log-scale.
logbeta  = THETAmatrix(check,1);
logdelta = THETAmatrix(check,2);
logalpha = THETAmatrix(check,3);
loggamma = THETAmatrix(check,4);
logtau = THETAmatrix(check,5);
logsigmabeta = THETAmatrix(check,6);
logsigmadelta = THETAmatrix(check,7);
logsigmaalpha = THETAmatrix(check,8);
logsigmaerror = THETAmatrix(check,9);

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

count_subj=1;
for subject = subjectsid'   % NEAT! if I take the transpose of subjectsid (so now it's a row vector) I can iterate through it
    sublogdata = alllogdata(alllogdata(:,3)==subject,:); %extract data pertaining the current subject
    n = size(sublogdata,1);  % number of measurements for the given subject
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


   lognoisyobs = log(allvolume) + sigmaerror*randn(size(allvolume,1),numsim);

   simsummaries = tumor_summaries(lognoisyobs,alllogdata(:,1),alllogdata(:,3));  % simulated summaries
   mean_simsummaries = mean(simsummaries,2);
   cov_simsummaries = cov(simsummaries');

   postpredcheck(check,:)  = mvnrnd(mean_simsummaries,cov_simsummaries);


end
save('postpredcheck','postpredcheck')

end
