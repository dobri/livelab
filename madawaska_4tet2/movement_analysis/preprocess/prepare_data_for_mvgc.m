function DATA = prepare_data_for_mvgc(DATA, dataTraj, headMark, plotting_flag)

% This function contains several preprocessing steps for the madawaska data
% The input DATA is all the ensemble data for the madawaska quartet
% The input dataTraj is the data you would like to process (e.g. X,
% X_detrended, the_trend, V, etc.
% The input headMark is a cell list of all the markers we want to analyze
% If you want to plot the downsampling, set plotting_flag=1. If not, set it to 0.

% First, I extract the head markers we are interested in (i.e. the four 
% head markers from each musician in the quartet)

% Then, I downsample the data to approximately 8hz

% Then, I convert all the data to z-scores

% Then, I take the average of the 4 head markers for each musician so I get
% one value for each musician

% written by Emily W, April 2020

%% Function starts here

% loop through all trials
for triali=1:length(DATA)
    
    % take out the data
    sf=DATA{triali}.sf; %sampling rate
    %recordedLength = (length(DATA{triali}.X))/sf; % length of data (in seconds)
    data = DATA{triali}.(dataTraj); % take out the actual data

    markLabels = DATA{triali}.col_names;

    % select only the head markers we want
    markInd=zeros(1,length(headMark));
    for n = 1:length(headMark) % find the 16 head markers
        markInd(n) = find(strcmp(markLabels, headMark{n}));
    end

    data=data(:,:,markInd);               

    % downsample the data
    timePoints=(1:length(data))'./sf;
    sf_target=8;
    ds_factor=round(sf/sf_target);
    length_data=size(data,1);

    step=(1:ds_factor:length_data)';
    windows_indices=[step min(step+ds_factor,length(data))];

    data_ds=zeros(size(windows_indices,1),3,16);
    timePoints_ds=zeros(size(windows_indices,1),1);

    for w=1:size(windows_indices,1)
        data_ds(w,:,:)=mean(data(windows_indices(w,1):windows_indices(w,2),:,:),1);
        timePoints_ds(w)=mean(timePoints(windows_indices(w,1):windows_indices(w,2)));
    end

    % optional plotting for verifying the downsampling
    % take just one dimension
    if plotting_flag==1
        data_ds_v=data_ds(:,1,1);
        data_v=data(:,1,1);

         sr=1/mean(diff(timePoints_ds)); % Your actual sr at the end.
         plot(timePoints,data_v)
         hold on
         plot(timePoints_ds,data_ds_v,'-or')
         hold off
    end
    
    % convert the data to z-scores
    if any(isnan(data_ds(:)))
    	mu=nanmean(data_ds);
    	sigma=nanstd(data_ds);
    	data_z=(data_ds-repmat(mu,length(data_ds),1))./repmat(sigma,length(data_ds),1);
    else
        data_z = zscore(data_ds);
    end              

    % take average of the 4 head markers for each musician
    data_cello = mean(data_z(:,:,1:4),3);
    data_viola = mean(data_z(:,:,5:8),3);
    data_violin1 = mean(data_z(:,:,9:12),3);
    data_violin2 = mean(data_z(:,:,13:16),3);
    
    data_avg=cat(3,data_cello,data_viola,data_violin1,data_violin2);

    % now take out just the anterior-posterior body sway
    data_AP=data_avg(:,2,:);
    
    % put processed data in a new field
    label=[dataTraj,'_processed'];
    DATA{triali}.(label)=data_AP;
end

