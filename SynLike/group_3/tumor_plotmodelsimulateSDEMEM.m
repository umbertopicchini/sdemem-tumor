function out = tumor_plotmodelsimulateSDEMEM(startstate,starttime,bigtheta,alllogdata,maxgrouptime)

% produces plots of simulated observations and removes points having
% tumor-logvolumes > log(1000) mm^3 as from experimental protocol

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


subjectsid = unique(alllogdata(:,3));
count_subj = 1;
for subject = subjectsid'   % NEAT! if I take the transpose of subjectsid (so now it's a row vector) I can iterate through it
    sublogdata = alllogdata(alllogdata(:,3)==subject,:); %extract data pertaining the current subject
    n = size(sublogdata,1);  % number of measurements for the given subject
   % v0 = sublogdata(1,2); 
    v0 = startstate(count_subj);
    count_subj = count_subj+1;
    times = sublogdata(:,1);
    stepsizes = diff(times);
    alpharand = mytruncgaussrandraw([0,1,alpha,sigmaalpha]);  %Umberto's custom made file
    betarand = beta + sigmabeta*randn;
    deltarand = delta + sigmadelta*randn;
    v0killed = alpharand*v0;
    v0surv = (1-alpharand)*v0;
    survived = zeros(n,1);  % n+1 rows because when plotting we also want to include the initial state 
    killed = survived;
    randn_std1 = randn(n,1);
    randn_std2 = randn(n,1);
    W1 = zeros(n,1);             % preallocate arrays ...
    W2 = W1;            
    W1(1) = sqrt(times(1)-starttime).*randn_std1(1);      
    W2(1) = sqrt(times(1)-starttime).*randn_std2(1);   
    W1(2:n) = sqrt(stepsizes).*randn_std1(2:n); 
    W2(2:n) = sqrt(stepsizes).*randn_std2(2:n);  
    volume = zeros(n,1);
     
     for ii=1:n
          if ii==1
              % the explicit solution of the state equations
              survived(1) = v0surv.*exp(betarand.*(times(1)-starttime)+gamma*W1(1));
              killed(1) = v0killed.*exp(-deltarand.*(times(1)-starttime)+tau*W2(1));
              volume(1) = survived(1) + killed(1);
          else
              % propagate forward from the particles having sampled indeces 
              survived(ii) = survived(ii-1).*exp(betarand.*stepsizes(ii-1)+gamma*W1(ii));
              killed(ii) = killed(ii-1).*exp(-deltarand.*stepsizes(ii-1)+tau*W2(ii));
              volume(ii) = survived(ii) + killed(ii);
          end
     end

    % and now add the initial assumed state (only for plotting purposes)
    volume = [v0;volume];

    % transform to log scale for consistency with data then add measurement
    % error
    individuallognoisyobs = log(volume) + sigmaerror*randn(n+1,1);  
    % now remove observations > log(1000) as from experimental protocol
    id_delete = find(individuallognoisyobs(2:end)>log(1000));
    if isempty(id_delete)
        id_keep = [1:length(individuallognoisyobs(2:end))]; % all observations are admissible
    else
        id_keep = [1:id_delete(1)]; % keep everything until the first one exceeding log(1000)
    end
    plot([starttime;times(id_keep)]*maxgrouptime,[individuallognoisyobs(1);individuallognoisyobs(id_keep+1)],'k-')  % plot individual simulated data vs unnormalised times
    hold on
%    plot(times*overallmaxtime,sublogdata(:,2),'bo-')
end
xlabel('days')
ylabel('log volume (mm^3)')
hold off
end

