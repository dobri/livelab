%% Where's my data?
plotting_video_flag = 0;
plotting_flag = 1;
data_folder = input('Type the path to your .tsv data files, i.e. ''D:\\str4tet_madawaska\\Data\'' or ''/Users/emilywood/Desktop/MATLAB/madawaska_4tet/data''\n>> ');
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
% fill the nan parts or cut them out, or the whole marker?cDATA{tr}.X

% Emily's matrix of NAN counts

% make a matrix
% the rows are the 83 markers
% the columns are the 12 trials
NAN_matrix = zeros(numel(DATA{1}.col_names), numel(filenames)); % This will store absolute count of NAN's
NAN_matrix_prop = NAN_matrix; % This will store the proportion of NAN's

for tr = 1:numel(DATA)
    for marker=1:size(DATA{tr}.X,3)
        NAN_matrix(marker, tr)=sum(isnan(DATA{tr}.X(:,1,marker))); % this is absolute count
        NAN_matrix_prop(marker,tr)=sum(isnan(DATA{tr}.X(:,1,marker))/32738); %proportion
    end
end

% ok there is a little bug in the code.. 
% the markers are different for files 1-8 files than for files 9-12 (corresponds to columns) 

% this is because the first 8 files are the full ensemble and the last 4 
% files are the soloists. There are 83 markers for the full ensemble but 
% less for each soloist (e.g. 19)

% e.g. compare DATA{1}.col_names to DATA{12}.col_names

% I don't think it matters though at this point. The soloist's markers will
% start from row 1, but once they run out, the rest of the values will be
% 0.

%% What other movement variables with reduced dimension should we think of?
% PCA? MSD? The thing from quaternions?
