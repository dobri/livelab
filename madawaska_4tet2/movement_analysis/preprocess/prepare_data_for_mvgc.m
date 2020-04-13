%% prepare_data_for_mvgc.m
% This script contains several preprocessing steps for the madawaska data

% First, I extract the head markers we are interested in (i.e. the four 
% head markers from each musician in the quartet)

% Then, I downsample the data to approximately 8hz

% Then, I convert all the data to z-scores

% Then, I take the average of the 4 head markers for each musician so I get
% one value for each musician

% Last, I create a matrix 'M' of desired data to take to the MVGC toolbox (i.e. 
% the anterior-posterior (AP) body sway of each musician's head markers in 
% the form: number_of_variables x number_of_observations x
% number_of_trials)

% written by Emily W, April 2020

%% Script starts here
% load DATA that Dobri sent me (from the MASTER_preprocess.m script), but
% only the 8 ensemble trials (cut out the solo data)

%load(DATA.mat)

% specify the names of the markers I want!
headMark = {'cellohat0','cellohat1','cellohat2','cellohat3','violahat0','violahat1','violahat2','violahat3',...
    'violin1hat0','violin1hat1','violin1hat2','violin1hat3','violin2hat0','violin2hat1','violin2hat2','violin2hat3'};

% now loop through all 8 trials
for triali=1:length(DATA)
    
    % take out the data
    sf=DATA{triali}.sf;
    recordedLength = (length(DATA{triali}.X))/sf; % length of data (in seconds)
    dataTraj = DATA{triali}.X; % take out the data
    markLabels = DATA{triali}.col_names;

    % select only the head markers
    markInd=zeros(1,length(headMark));
    for n = 1:length(headMark) % find the 16 head markers
        markInd(n) = find(strcmp(markLabels, headMark{n}));
    end

    dataTraj=dataTraj(:,:,markInd);

    % downsample the data
    timePoints=(1:length(dataTraj))'./sf;
    sf_target=8;
    ds_factor=round(sf/sf_target);
    length_data=size(dataTraj,1);

    step=(1:ds_factor:length_data)';
    windows_indices=[step min(step+ds_factor,length(dataTraj))];

    dataTraj_ds=zeros(size(windows_indices,1),3,16);
    timePoints_ds=zeros(size(windows_indices,1),1);

    for w=1:size(windows_indices,1)
        dataTraj_ds(w,:,:)=mean(dataTraj(windows_indices(w,1):windows_indices(w,2),:,:),1);
        timePoints_ds(w)=mean(timePoints(windows_indices(w,1):windows_indices(w,2)));
    end
    
    % optional plotting for verifying the downsampling
    % take just one dimension
    dataTraj_ds_v=dataTraj_ds(:,1,1);
    dataTraj_v=dataTraj(:,1,1);

    sr=1/mean(diff(timePoints_ds)); % Your actual sr at the end.
    plot(timePoints,dataTraj_v)
    hold on
    plot(timePoints_ds,dataTraj_ds_v,'-or')
    hold off
    
    % convert the data to z-scores
    dataTraj_z = zscore(dataTraj_ds);
    DATA{triali}.Z=dataTraj_z;

    % take average of the 4 head markers for each musician
    data_cello = mean(dataTraj_z(:,:,1:4),3);
    data_viola = mean(dataTraj_z(:,:,5:8),3);
    data_violin1 = mean(dataTraj_z(:,:,9:12),3);
    data_violin2 = mean(dataTraj_z(:,:,13:16),3);
    
    dataTraj_avg=cat(3,data_cello,data_viola,data_violin1,data_violin2);
    DATA{triali}.AVG=dataTraj_avg;

    % now take out just the anterior-posterior body sway
    dataTraj_AP=dataTraj_avg(:,2,:);
    DATA{triali}.AP=dataTraj_AP;

end

%% Create a matrix for MVGC

% now I need a matrix to take over to MVGC toolbox!

% first, let's see how long each trial is
trialLengths=zeros(1,length(DATA));
for i=1:length(DATA)
    trialLengths(1,i)=length(DATA{i}.AP);
end

minVal=min(trialLengths);
% ok I'm going to cap all of them at minVal (2513) because that's the shortest

% now create the matrix M for MVGC
nvars=size(DATA{1}.AP,3);
nobs=minVal;
ntrials=length(DATA);

M=zeros(nvars,nobs,ntrials);

for tri=1:ntrials
    AP=permute(DATA{tri}.AP,[3,1,2]);
    M(:,:,tri)=AP(:,1:minVal);
end

%% Plot the matrix to make sure it looks ok

for i=1:size(M,3) %loop through all 8 trials
    for j=1:size(M,1) % loop through all 4 musicians
        plot(M(j,:,i))
        pause
    end
end

%% Save my matrix and updated DATA
save('M.mat','M')
save('DATA.mat','DATA')



