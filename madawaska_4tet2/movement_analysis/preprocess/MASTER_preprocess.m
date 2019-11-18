%% Where's my data?
plotting_video_flag = 0;
plotting_flag = 0;
data_folder = input('Type the path to your .tsv data files, i.e. ''D:\\str4tet_madawaska\\Data\''>>' );
fprintf('%s\n',['We''ll be looking inside ' data_folder])
filenames = dir(fullfile(data_folder,'*.tsv'));


%% Import the QTM data saved as text files.
for tr = 1:numel(filenames)
    filename = filenames(tr); %filename = 'piece1_solo1.tsv'
    fprintf('%s\n',filename.name)
    DATA{tr} = import_tsv_from_qtm_to_matlab(filename);
end


%% To confirm the data looks OK set the flag to 1 and run an animated plot.
for tr = 1:numel(DATA)
    fprintf('%s\n',filenames(tr).name)
    if plotting_video_flag == 1
        plot_animated_in_3d(DATA{tr})
    end
end


%% Here we do cleaning, nan-filling, detrending, zero-center, etc..
% The time series should look reasonably stationary from here on.
% There's the option to plot each marker's raw and detrended time series.
for tr = 1:numel(DATA)
    % De-trend, de-nan, de-artefact.
    % Fill gaps by linearly interpolating, good enough for small gaps.
    DATA{tr}.X_detrended = zeros(size(DATA{tr}.X,1),size(DATA{tr}.X,2),size(DATA{tr}.X,3));
    DATA{tr}.The_trend = nan(size(DATA{tr}.X,1),size(DATA{tr}.X,2),size(DATA{tr}.X,3));
    for marker=1:size(DATA{tr}.X,3)
        [DATA{tr}.X_detrended(:,:,marker), DATA{tr}.The_trend(:,:,marker)] = remove_nonstation_and_nans_with_splines(DATA{tr}.X(:,:,marker),DATA{tr}.sf,(1:size(DATA{tr}.X,1))'./DATA{tr}.sf,[],[],plotting_flag);
        if plotting_flag == 1
            text(.45,-.3,[DATA{tr}.filename(~(double(DATA{tr}.filename)==95)) ' - ' DATA{tr}.col_names{marker}],'units','normalized')
            pause
        end
    end
end


%% Get proportion of NaNs for each participant and marker.
% What to do with nan gaps? Decide depending on how broken is the data, 
% fill the nan parts or cut them out, or the whole marker?


%% What other movement variables with reduced dimension should we think of?
% PCA? MSD? The thing from quaternions?