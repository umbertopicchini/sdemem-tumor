

rng(100)


problem = 'tumor';
group = 3; % select the treatment group to fit 

%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
allRAWdata = dlmread('tumorlong_nofollowup.txt',' ',1,0); % read data having follow ups manually removed
alllogdata = zeros(size(allRAWdata,1),size(allRAWdata,2));
% let's rearrange the columns the way we prefer
alllogdata(:,1) = allRAWdata(:,1); % times
alllogdata(:,2) = log(allRAWdata(:,3)); % volumes in log-scale, store them in column 2
alllogdata(:,3) = allRAWdata(:,2);  % the subject id, goes in column 3
alllogdata(:,4) = allRAWdata(:,4);  % the experimental group
overallmaxtime = max(alllogdata(:,1)); % the maximum observed time across ALL groups

% extract the desired group from the dataset
logsubdata = alllogdata(alllogdata(:,4)==group,1:3);  % the dataset for subjects in group "group"

% remove dying mouse from group 1 and subjects with too few observations
if group==1
    logsubdata = logsubdata(ismember(logsubdata(:,3), [11,13,15,17,18] ) ,:);
end

%maxgrouptime = max(logsubdata(:,1));   % max time for the given group; only to be used when plotting
logsubdata(:,1) = logsubdata(:,1)/overallmaxtime; % scale all times by the maximum observed time (for numerical stability)
%plotdata(alllogdata,group); % let's look at the log-data for the chosen group 

numsubjects = length(unique(logsubdata(:,3)));
startstate = zeros(1,numsubjects);
subjectsid = unique(logsubdata(:,3));
count_subj = 1;
for subject = subjectsid'   % NEAT! if I take the transpose of subjectsid (so now it's a row vector) I can iterate through it
    subjectlogdata = logsubdata(logsubdata(:,3)==subject,:);
    startstate(count_subj) = exp(subjectlogdata(1,2)); % [mm^3] hypothetical subject's starting volume (NOT LOG VOLUME). Only useful to initialize the simulation
    count_subj = count_subj+1;
end

% if it is an experimental group (i.e. if group < 5) remove data taken at
% day 1 and day 4 (because these are days when mice get injections --> it
% takes time to observe an effect). Hence when "group<5" the resulting state-space model
% only consider measurement times from day 6 onward
if group < 5
   logsubdata = logsubdata(logsubdata(:,1)>=6/overallmaxtime,:);
end

starttime = 0; %i.e. day 0


%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


%:: parameters starting values
%         log(beta),   logdelta  logalpha   loggamma    logtau,   logsigmabeta,  logsigmadelta   logsigmaalpha   logsigmaerror
parbase = [1.6          1.6       -0.36        0         0       -0.7              -0.7               -2.3                0    ];
parmask = [1             1           1         1         1         1                 1                  1                 1   ];  % 1 for parameters to be estimated, 0 otherwise
bigtheta = parbase;


R_mcmc = 20000; % number of MCMC iterations
step_rw = [0.1 0.1 0.01 0.1 0.1 0.1 0.1 0.1 0.1];%initial values for the standard deviations of adaptive Gaussian random walk 
covariance_refresh = 100; % the frequency of update for the covariance matrix in adaptive Metropolis 
numparticles = 2000; % number of particles for the auxiliary particle filter (denoted L in the paper)
numAPFsim = 5;   % number of particles required to compute the means in the first stage of the auxiliary particle filter (denoted L2 in the paper)


THETAmatrix = pmcmc_sdmem_tumor_APF(problem,group,logsubdata,startstate,starttime,bigtheta,parmask,parbase,R_mcmc,step_rw,covariance_refresh,numparticles,numAPFsim);

% posterior means reported in the paper are obtained with
mean(exp(THETAmatrix(10000:end,:)));


meantheta = mean(THETAmatrix(10000:end,:));
% plot simulated trajectories
tumor_plotmodelsimulateSDEMEM(startstate,starttime,meantheta,logsubdata,overallmaxtime);

