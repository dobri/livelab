function DATA = prepare_data_for_mvgc(DATA, headMark, plotting_flag)

% This function contains several preprocessing steps for the madawaska data
% The input DATA is all the ensemble data for the madawaska quartet
% The input headMark is a cell list of all the markers we want to analyze
% If you want to plot everything, set plotting_flag=1. If not, set it to 0.

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

%% Function starts here

% loop through all trials
for triali=1:length(DATA)
    
    % take out the data
    sf=DATA{triali}.sf; %sampling rate
    %recordedLength = (length(DATA{triali}.X))/sf; % length of data (in seconds)
    dataTraj = DATA{triali}.X; % take out the actual data
    markLabels = DATA{triali}.col_names;

    % select only the head markers we want
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
    if plotting_flag==1
        dataTraj_ds_v=dataTraj_ds(:,1,1);
        dataTraj_v=dataTraj(:,1,1);

         sr=1/mean(diff(timePoints_ds)); % Your actual sr at the end.
         plot(timePoints,dataTraj_v)
         hold on
         plot(timePoints_ds,dataTraj_ds_v,'-or')
         hold off
    end
    
    % convert the data to z-scores
    if any(isnan(dataTraj_ds(:)))
    	mu=nanmean(dataTraj_ds);
    	sigma=nanstd(dataTraj_ds);
    	dataTraj_z=(dataTraj_ds-repmat(mu,length(dataTraj_ds),1))./repmat(sigma,length(dataTraj_ds),1);
    else
        dataTraj_z = zscore(dataTraj_ds);
    end              

    % take average of the 4 head markers for each musician
    data_cello = mean(dataTraj_z(:,:,1:4),3);
    data_viola = mean(dataTraj_z(:,:,5:8),3);
    data_violin1 = mean(dataTraj_z(:,:,9:12),3);
    data_violin2 = mean(dataTraj_z(:,:,13:16),3);
    
    dataTraj_avg=cat(3,data_cello,data_viola,data_violin1,data_violin2);

    % now take out just the anterior-posterior body sway
    dataTraj_AP=dataTraj_avg(:,2,:);
    DATA{triali}.AP=dataTraj_AP;

end

