function allsummaries = tumor_summaries(logvolume,timesvec,subjectsvec)
% Summary statistics for Bayesian synthetic likelihood calculations
%
%Input: - log-volumes (observed or simulated)
%       - sampling times
%       - ID of the subjects
% Output: a vector containing 5*numsubjects+3 statistics, i.e. observed or simulated synthetic summaries, corresponding to 
%                         5 individual summaries for each subject, and additional 3 summaries that are common/shared
%                         between all subjects. E.g. group 3 has 5 subjects, therefore it outputs 28 values.

numsim = size(logvolume,2);
summaries = zeros(5,numsim); 
allsummaries = [];
subjectsid = unique(subjectsvec);

% summaries computed on each trajectory separately (i.e. individual
% summaries)
for subject = subjectsid'   % NEAT! if I take the transpose of subjectsid (so now it's a row vector) I can iterate through it
    sublogvolume = logvolume(subjectsvec==subject,:); %extract data pertaining the current subject
    subtimes = timesvec(subjectsvec==subject,:);
    summaries(1,:) = mad(sublogvolume,0,1);
    summaries(2,:) = (sublogvolume(end,:)-sublogvolume(1,:))./(subtimes(end,:)-subtimes(1,:));
    summaries(3,:) = sublogvolume(1,:);
    summaries(4,:) = sublogvolume(2,:); 
    % least squares slopes
    summaries(5,:) = sum( (sublogvolume(1:end-1,:) - repmat(mean(sublogvolume(1:end-1,:),1),size(sublogvolume,1)-1,1) ) .* (sublogvolume(2:end,:) - repmat(mean(sublogvolume(2:end,:),1),size(sublogvolume,1)-1,1) ) ,1) ./ sum((sublogvolume(1:end-1,:) - repmat(mean(sublogvolume(1:end-1,:),1),size(sublogvolume,1)-1,1) ).^2, 1); 
    allsummaries = [allsummaries;summaries];
end

% augment with summaries computed on the whole "ensemble" 
allsummaries = [allsummaries; (mad(logvolume(timesvec==timesvec(1),:),0,1))]; % computes MAD on the first time point between all trajectories
allsummaries = [allsummaries; (mad(logvolume(timesvec==timesvec(2),:),0,1))]; % computes MAD on the third time point between all trajectories
allsummaries = [allsummaries; (mad(logvolume(timesvec==timesvec(end),:),0,1))];
