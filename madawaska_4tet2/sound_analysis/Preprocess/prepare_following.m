function data_z = prepare_following(X, ds_target, sf, plotting_flag)

%% This function contains the preprocessing steps for the musical following study
% Downsample a time series and convert it to z-scores
% The input is a raw time series as a (nobs x 1) matrix
% To plot the downsampling, set plotting_flag = 1. If not, set to 0.



%% DOWNSAMPLE
timePoints = (1:length(X))'./sf; % # of seconds in data
ds_factor = round(sf/ds_target); % 
length_data = size(X,1); % # of samples in data

% Make a matrix with two columns that indexes the beginning (col. 1
% and the end (col. 2) of each window
step = (1:ds_factor:length_data)'; % Array with steps of size (ds_factor)
windows_indices = [step min(step+ds_factor,length(X))];

data_ds=zeros(size(windows_indices,1),1);
timePoints_ds=zeros(size(windows_indices,1),1);

for w=1:size(windows_indices,1)
    % For each time point, take the mean of the two time points
    % that correspond to the start and end of the window
    data_ds(w,:,:)=mean(X(windows_indices(w,1):windows_indices(w,2),:,:),1);
    
    % Index the timePoints array using the windows_indices array
    timePoints_ds(w)=mean(timePoints(windows_indices(w,1):windows_indices(w,2)));
end


% Optional plotting to check downsampling
if plotting_flag == 1
    % data_ds_v = data_ds(:,1)
    %sr = 1/mean(diff(timePoints_ds));
    plot(timePoints_ds,data_ds,'or')
    hold off
end



%% Z-SCORES
% Convert the data to z-scores
if any(isnan(data_ds(:)))
    mu=nanmean(data_ds);
    sigma=nanstd(data_ds);
    data_z=(data_ds-repmat(mu,length(data_ds),1))./repmat(sigma,length(data_ds),1);
else
    data_z = zscore(data_ds);
end





