function M = create_matrix_following(DATA, plotting_flag)

%% Create a matrix for the MVGC toolbox
% This function takes data from the amplitude envelopes of two time series
% and puts it in the correct form for the MVGC toolbox.
% This data is comprised of two variables: a musician's performance and the
% recording they followed along with. The recording is the same for each trial. 

%% So for each trial,
% var #1 is the performance and var #2 is the recording.
% nobs is the length of that trial. 



%% Where's my data?
trialLengths = zeros(1,length(DATA)); % How long is each trial?

for triali = 1:length(DATA) % number of trials

    trialLengths(1,triali) = length(DATA{triali}.ds_100);
    
end

minVal = min(trialLengths); % Find the shortest trial

nvars = size(DATA{1}.ds_100,3); % We can just use trial 1 because all trials
% should be the same length (?)
nobs = minVal;
ntrials = length(DATA);

M = zeros(nvars,nobs,ntrials);

for tri=1:ntrials
    data = permute(DATA{tri}.ds_100,[3,1,2]); %(?)
    M(:,:,tri) = data(:,1:minVal);
end



%% Optional - plot the matrix to make sure it looks ok 
if plotting_flag ==1

    for triali=1:size(M,3) %loop through all trials
        for j=1:size(M,1) % loop through both vars (recording and performance)
            plot(M(j,:,triali))
            disp([triali, j])
            pause
        end
    end
end
