function DATA = prepare_data_for_mvgc_simpler(DATA, dataTraj, dim, plotting_flag)
% This function contains several preprocessing steps for the madawaska data
% The input DATA is all the ensemble data for the madawaska quartet
% The input dataTraj is the data you would like to process (e.g. X,
% X_detrended, the_trend, V, etc.
% Then, I downsample the data to approximately 8hz
% Then, I convert all the data to z-scores
% written by Emily W, April 2020

%% Function starts here

% loop through all trials
for triali=1:length(DATA)
    % take out the data
    sf=DATA{triali}.sf; %sampling rate
    if ndims(DATA{triali}.(dataTraj))==3
        data = squeeze(DATA{triali}.(dataTraj)(:,dim,:)); % take out the actual data
    else
        data = DATA{triali}.(dataTraj); % take out the actual data
    end
    
    % downsample the data
    timePoints=(1:length(data))'./sf;
    sf_target=8;
    ds_factor=round(sf/sf_target);
    length_data=size(data,1);
    
    step=(1:ds_factor:length_data)';
    windows_indices=[step min(step+ds_factor,length(data))];
    
    data_ds=zeros(size(windows_indices,1),size(data,2));
    timePoints_ds=zeros(size(windows_indices,1),1);
    
    for w=1:size(windows_indices,1)
        data_ds(w,:)=mean(data(windows_indices(w,1):windows_indices(w,2),:),1);
        timePoints_ds(w)=mean(timePoints(windows_indices(w,1):windows_indices(w,2)));
    end
    
    % optional plotting for verifying the downsampling
    % take just one dimension
    if plotting_flag==1
        plot(timePoints,data)
        hold on
        plot(timePoints_ds,data_ds,'o')
        hold off
        pause
    end
    
    % convert the data to z-scores
    if any(isnan(data_ds(:)))
        mu=nanmean(data_ds);
        sigma=nanstd(data_ds);
        data_z=(data_ds-repmat(mu,length(data_ds),1))./repmat(sigma,length(data_ds),1);
    else
        data_z = zscore(data_ds);
    end
    
    % put processed data in a new field
    label = [dataTraj,'_processed'];
    DATA{triali}.(label) = data_z;
end
