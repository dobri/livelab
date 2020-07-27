clear


%% Where's my data?
if strcmp(computer,'GLNXA64')
    data_folder = input('Where is my fokin money? '); % https://www.youtube.com/watch?v=pslUqRMPthk
else
    data_folder = '/Users/emilywood/Desktop/MATLAB/trainorlab/madawaska_4tet/data/ensemble';
end
which_piece = ['piece' num2str(input('Which of the two pieces? Type a number: '),'%1.0f')];
plotting_flag = input('To visualize data along the way? Type 1 for yes: ');


%% Import the QTM data saved as text files.
% First, import all the 4 solo recordings as trial 1 (surrogate coupling).
filenames = dir(fullfile(data_folder,[which_piece '_sol*.tsv']));
for tr = 1:numel(filenames)
    filename = filenames(tr);
    fprintf('%s\n',filename.name)
    % Why? This works on a Mac?
    if ~strcmp(computer,'GLNXA64')
        filename=filename.name;
    end
    DATA0{tr} = import_tsv_from_qtm_to_matlab(filename);
end

shortest_trial_len = inf;
for tr = 1:numel(filenames)
    shortest_trial_len = min(shortest_trial_len,size(DATA0{tr}.X,1));
end

clear DATA
tr=1;
DATA{1}.filename = DATA0{tr}.filename(1:end-5);
DATA{1}.col_names = DATA0{tr}.col_names;
DATA{1}.X = DATA0{tr}.X(1:shortest_trial_len,:,:);
DATA{1}.sf = DATA0{tr}.sf;
DATA{1}.ensemble_condition = 0;
for tr = 2:numel(filenames)
    DATA{1}.col_names = [DATA{1}.col_names DATA0{tr}.col_names];
    DATA{1}.X(:,:,(end+1):(end+size(DATA0{tr}.X,3))) = DATA0{tr}.X(1:shortest_trial_len,:,:);
end


% Second, import the ensemble recordings as trials 2+.
filenames = dir(fullfile(data_folder,[which_piece '_ens*.tsv']));
for tr = 1:numel(filenames)
    filename = filenames(tr); %filename = 'piece1_solo1.tsv'
    fprintf('%s\n',filename.name)
    if ~strcmp(computer,'GLNXA64')
        filename=filename.name;
    end
    DATA{tr+1} = import_tsv_from_qtm_to_matlab(filename);
    DATA{tr+1}.ensemble_condition = 1;
end
bodies_labels = {'violin1','violin2','viola','cello'};


%% Should we rotate the bodies so that the AP and ML axes correspond to Y 
% and X axes in the data frame, respectively?
target_vector = [-1 0];
dims = [1 2];
reference_markers = {'backl','backr'};
for tr = 1:numel(DATA)
    fprintf('%s\n',DATA{tr}.filename)
    DATA{tr} = rotate_to_target_vector(DATA{tr},target_vector,dims,bodies_labels,reference_markers,plotting_flag);
    if plotting_flag==1;pause;end
end


%% Get proportion of NaNs for each participant and marker.
% Trials are cells, markers are columns.
for tr = 1:numel(DATA)
    for marker=1:size(DATA{tr}.X,3)
        %DATA{tr}.missing_count(1,marker) = sum(isnan(DATA{tr}.X(:,1,marker))); % this is absolute count
        DATA{tr}.missing_prop(1,marker) = sum(isnan(DATA{tr}.X(:,1,marker)))/size(DATA{tr}.X(:,1,marker),1); %proportion
    end
end
% What to do with nan gaps? Decide depending on how broken is the data. 
% Fill the nan parts or cut them out, or the whole marker?


%% Show the amount of missing data using figures and printing to screen.
% Let's first focus on head sway as a proxy for upper body gross movement.
wanted_head_markers = {'hat0','hat1','hat2','hat3'};
if plotting_flag == 1
    for tr = 1:numel(DATA)
        for body = 1:4
            for marker = 1:numel(wanted_head_markers)
                for m = 1:numel(DATA{tr}.col_names)
                    if strcmp(DATA{tr}.col_names{m},[bodies_labels{body} wanted_head_markers{marker}])
                        plot(body*10+marker,DATA{tr}.missing_prop(m),'s');
                        fprintf('%6.4f,',DATA{tr}.missing_prop(m))
                        text(body*10+marker,.1,DATA{tr}.col_names{m},'Rotation',90)
                        hold on
                    end
                end
            end
        end
        fprintf('\n')
        set(gca,'XTick',[])
        hold off
        ylim([0 1])
        pause
    end
end
% It seems that head markers were perfectly recorded, at least for piece1.
% The head markers have some missing values for trials 6 and 8 in piece 2.


%% To confirm the data looks OK set the flag to 1 and run an animated plot.
plotting_video_flag = 1;
if plotting_video_flag == 1
    for tr = 1:numel(DATA)
        fprintf('%s\n',DATA{tr}.filename)
        plot_animated_in_3d(DATA{tr}.X)
        pause
    end
end


%% Nicely looking 3D multi-trial, multi-body video.
wanted_head_markers_all = {'hat0','scroll1','backl','spine1','backr','elbowl','wristl','handl','elbowr','wristr','handr','bow1','bow2','boutr','boutl'};
% plot_4x8_animated_in_3d(DATA,bodies_labels,wanted_head_markers_all,100);


%% Here we do cleaning, nan-filling, detrending, zero-center, etc..
% The time series should look reasonably stationary from here on.
% There's the option to plot each marker's raw and detrended time series.
for tr = 1:numel(DATA)
    % De-trend, de-nan, de-artefact.
     %Fill gaps by linearly interpolating, good enough for small gaps.
     DATA{tr}.X_detrended = zeros(size(DATA{tr}.X,1),size(DATA{tr}.X,2),size(DATA{tr}.X,3));
     DATA{tr}.the_Trend = nan(size(DATA{tr}.X,1),size(DATA{tr}.X,2),size(DATA{tr}.X,3));
     for marker=1:size(DATA{tr}.X,3)
        [DATA{tr}.X_detrended(:,:,marker), DATA{tr}.the_Trend(:,:,marker)] = remove_nonstation_and_nans_with_splines(DATA{tr}.X(:,:,marker),DATA{tr}.sf,(1:size(DATA{tr}.X,1))'./DATA{tr}.sf,[],[],plotting_flag);
        if plotting_flag == 1
            text(.45,-.3,[DATA{tr}.filename(~(double(DATA{tr}.filename)==95)) ' - ' DATA{tr}.col_names{marker}],'units','normalized')
            pause
        end
    end
end


%% Speed of head movement. MSD?
for tr = 1:numel(DATA)
    for b = 1:numel(bodies_labels)
        % Find the indices of the needed markers.
        marker_index = zeros(1,numel(wanted_head_markers));
        for m = 1:numel(DATA{tr}.col_names)
            for marker = 1:numel(marker_index)
                if strcmp(DATA{tr}.col_names{m},[bodies_labels{b} wanted_head_markers{marker}])
                    marker_index(marker) = m;
                end
            end
        end
        
        % If no markers for this body are found.
        if all(marker_index==0)
            continue
        end
        
        % Extract the 3D values of these markers.
        X = DATA{tr}.X(:,:,marker_index);
        
        % We don't actually need all four head markers. %X = X(:,:,1);
        % If there are issues with missing, you can try averaging across
        % them to smooth the data a little.
        X = nanmean(X,3);
        
        % Also, verify that some markers do not disappear a lot, which 
        % would make the average jump up/down by a few mm.
        % Some glitches in the raw 3D cause huge differences in v.
        % It's easier to pre-smooth X before v, rather than to clean v.
        v_temp = [0;dot(diff(X),diff(X),2).^.5./(1/DATA{tr}.sf)];
        X(logical((v_temp>3e2).*([0;abs(diff(v_temp))>1e2])),:)=nan;
        for d=1:3
            X(:,d) = fill_nans_by_lin_interp(X(:,d));
            % Smooth by a third of a second. With the mov ave method this
            % kills everything above 3 Hz. With sgolay above 5 Hz.
            X(:,d) = smooth(X(:,d),round(DATA{tr}.sf/3),'sgolay');
        end
        
        % The speed, excluding some in and out window.
        v = [0;dot(diff(X),diff(X),2).^.5./(1/DATA{tr}.sf)];
        v([1:300 end-299:end]) = 0;
        DATA{tr}.V(:,b) = v;
        
        % Verify that there aren't spikes and nans remaining by accident.
        % Did we do a good job cleaning and filtering without killing v?
        clf
        subplot(1,3,1)
        plot(DATA{tr}.X(:,:,marker_index(2)));
        subplot(1,3,2)
        plot(X);
        subplot(1,3,3)
        plot([v_temp DATA{tr}.V(:,b)])
        pause
    end
end

%% Save the whole DATA structure to a mat file.
% save(fullfile(data_folder,['DATA_piece' which_piece(end) '.mat']),'DATA','-v7.3')

%% START HERE IF YOU ALREADY HAVE DATA FILES FROM DOBRI
% first, we'll define some variables from above that we haven't run yet
bodies_labels = {'violin1','violin2','viola','cello'};
plotting_flag=1;

% load the first datafile
load('DATA_piece1.mat');

% rename it DATA1 (because loading DATA_piece2 will overwrite it)
DATA1=DATA;

% Take out first trial (because this is the solo trial. I'm working with 
% ensemble data first)
DATA1(1)=[];

% now load the data from piece 2
load('DATA_piece2.mat');

% rename it DATA2 
DATA2=DATA;

% Take out first trial for the same reason stated above
DATA2(1)=[];

%% Show the amount of missing data using figures and printing to screen.
%I already know that DATA1 was recording perfectly. So let's check DATA2
wanted_head_markers = {'hat0','hat1','hat2','hat3'};
if plotting_flag == 1
    for tr = 1:numel(DATA2)
        for body = 1:4
            for marker = 1:numel(wanted_head_markers)
                for m = 1:numel(DATA2{tr}.col_names)
                    if strcmp(DATA2{tr}.col_names{m},[bodies_labels{body} wanted_head_markers{marker}])
                        plot(body*10+marker,DATA2{tr}.missing_prop(m),'s');
                        fprintf('%6.4f,',DATA2{tr}.missing_prop(m))
                        text(body*10+marker,.1,DATA2{tr}.col_names{m},'Rotation',90)
                        hold on
                    end
                end
            end
        end
        fprintf('\n')
        set(gca,'XTick',[])
        hold off
        ylim([0 1])
        pause
    end
end
% The head markers have some missing values for trials 6 and 8 in piece 2.

%% FILL IN NANS in DATA2
%Fill gaps with linear interpolation in the X vector
for tr = 1:numel(DATA2)
    DATA2{tr}.X_old=DATA2{tr}.X;
    DATA2{tr}.X = zeros(size(DATA2{tr}.X,1),size(DATA2{tr}.X,2),size(DATA2{tr}.X,3));
     for marker=1:size(DATA2{tr}.X,3)
        DATA2{tr}.X(:,:,marker) = fill_nans_with_splines_X(DATA2{tr}.X_old(:,:,marker));
    end
end

%check to see if it worked
headMark = {'cellohat0','cellohat1','cellohat2','cellohat3','violahat0','violahat1','violahat2','violahat3',...
    'violin1hat0','violin1hat1','violin1hat2','violin1hat3','violin2hat0','violin2hat1','violin2hat2','violin2hat3'};

for tr = 1:numel(DATA2)
    fprintf('Trial %1.0f\n',tr)
    for marker=1:size(DATA2{tr}.X,3)
        if any(strcmp(headMark,DATA2{tr}.col_names{marker}))
            fprintf('%6.4f\n',sum(isnan(DATA2{tr}.X(:,1,marker)))/size(DATA2{tr}.X(:,1,marker),1)); %proportion of nans. should all be 0
        end
    end
end

% Redo cleaning, nan-filling, detrending, zero-center, etc. for DATA2
plotting_flag=0;
for tr = 1:numel(DATA2)
    % De-trend, de-nan, de-artefact.
     %Fill gaps by linearly interpolating, good enough for small gaps.
     DATA2{tr}.X_detrended = zeros(size(DATA2{tr}.X,1),size(DATA2{tr}.X,2),size(DATA2{tr}.X,3));
     DATA2{tr}.the_Trend = nan(size(DATA2{tr}.X,1),size(DATA2{tr}.X,2),size(DATA2{tr}.X,3));
     for marker=1:size(DATA2{tr}.X,3)
        [DATA2{tr}.X_detrended(:,:,marker), DATA2{tr}.the_Trend(:,:,marker)] = remove_nonstation_and_nans_with_splines(DATA2{tr}.X(:,:,marker),DATA2{tr}.sf,(1:size(DATA2{tr}.X,1))'./DATA2{tr}.sf,[],[],plotting_flag);
        if plotting_flag == 1
            text(.45,-.3,[DATA2{tr}.filename(~(double(DATA2{tr}.filename)==95)) ' - ' DATA2{tr}.col_names{marker}],'units','normalized')
            pause
        end
    end
end

% redo speed of head movement
for tr = 1:numel(DATA2)
    for b = 1:numel(bodies_labels)
        % Find the indices of the needed markers.
        marker_index = zeros(1,numel(wanted_head_markers));
        for m = 1:numel(DATA2{tr}.col_names)
            for marker = 1:numel(marker_index)
                if strcmp(DATA2{tr}.col_names{m},[bodies_labels{b} wanted_head_markers{marker}])
                    marker_index(marker) = m;
                end
            end
        end
        
        % If no markers for this body are found.
        if all(marker_index==0)
            continue
        end
        
        % Extract the 3D values of these markers.
        X = DATA2{tr}.X(:,:,marker_index);
        
        % We don't actually need all four head markers. %X = X(:,:,1);
        % If there are issues with missing, you can try averaging across
        % them to smooth the data a little.
        X = nanmean(X,3);
        
        % Also, verify that some markers do not disappear a lot, which 
        % would make the average jump up/down by a few mm.
        % Some glitches in the raw 3D cause huge differences in v.
        % It's easier to pre-smooth X before v, rather than to clean v.
        v_temp = [0;dot(diff(X),diff(X),2).^.5./(1/DATA2{tr}.sf)];
        X(logical((v_temp>3e2).*([0;abs(diff(v_temp))>1e2])),:)=nan;
        for d=1:3
            X(:,d) = fill_nans_by_lin_interp(X(:,d));
            % Smooth by a third of a second. With the mov ave method this
            % kills everything above 3 Hz. With sgolay above 5 Hz.
            X(:,d) = smooth(X(:,d),round(DATA2{tr}.sf/3),'sgolay');
        end
        
        % The speed, excluding some in and out window.
        v = [0;dot(diff(X),diff(X),2).^.5./(1/DATA2{tr}.sf)];
        v([1:300 end-299:end]) = 0;
        DATA2{tr}.V(:,b) = v;
        
        % Verify that there aren't spikes and nans remaining by accident.
        % Did we do a good job cleaning and filtering without killing v?
        if plotting_flag==1
            clf
            subplot(1,3,1)
            plot(DATA2{tr}.X(:,:,marker_index(2)));
            subplot(1,3,2)
            plot(X);
            subplot(1,3,3)
            plot([v_temp DATA2{tr}.V(:,b)])
            pause
        end
    end
end

%% Now prepare data for gc analysis 
% specify the names of the markers I want to analyze
headMark = {'cellohat0','cellohat1','cellohat2','cellohat3','violahat0','violahat1','violahat2','violahat3',...
    'violin1hat0','violin1hat1','violin1hat2','violin1hat3','violin2hat0','violin2hat1','violin2hat2','violin2hat3'};

% From our three data trajectories 'X', 'X_detrended' and 'the_Trend', take 
% the desired markers (head markers), downsample the data, convert to  
% z-scores, and take the average of the 4 head markers for each musician.

dataTrajs={'X','X_detrended','the_Trend'};
for i = 1:length(dataTrajs)
    DATA1=prepare_data_for_mvgc(DATA1,dataTrajs{i},headMark,plotting_flag);
    DATA2=prepare_data_for_mvgc(DATA2,dataTrajs{i},headMark,plotting_flag);
end

% Our DATA variables have new structs called X_processed, X_detrended_processed
% and the_Trend_processed, which correspond to the anterior-posterior body 
% sway of the 4 musicians in each of these trajectories.

% Get the data into a matrix form for the MVGC toolbox
dataTrajs={'X_processed','X_detrended_processed','the_Trend_processed'};
for i = 1:length(dataTrajs)
  D{1}.(dataTrajs{i})=create_matrix_for_mvgc(DATA1, dataTrajs{i}, plotting_flag);
  D{2}.(dataTrajs{i})=create_matrix_for_mvgc(DATA2, dataTrajs{i}, plotting_flag);
end

% Now we have a NEW data variable called D.
% There are TWO structs in it corresponding to piece 1 and piece 2
% WIthin these structs, we have a matrix for mvgc for each of our three
% data trajectories, including X, X_detrended and the_Trend

% Now take this variable D over to the mvgc toolbox for gc!
save('D','D')


%% ... What other movement variables with reduced dimension should we use?
% v AP sway (front-back). (Y-axis in the rotated data).
% v MSD or speed. (Do that!)
% x PCA? (Probably not necessary here to fine-tune beyond ML and AP. We know that PC1 and PC2 will correspond closely to the AP and ML axes.)
% x The thing from quaternions? (No time for that.)
% ACCELERATION!

