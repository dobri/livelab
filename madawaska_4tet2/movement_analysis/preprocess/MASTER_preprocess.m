clear

%% Where's my data?
plotting_flag = 0;
if strcmp(computer,'GLNXA64')
    data_folder = input('Where is my fokin money? ');
else
    data_folder = '/Users/emilywood/Desktop/MATLAB/trainorlab/madawaska_4tet/data/ensemble';
end
filenames = dir(fullfile(data_folder,'piece1*.tsv'));


%% Import the QTM data saved as text files.
for tr = 1:numel(filenames)
    filename = filenames(tr); %filename = 'piece1_solo1.tsv'
    fprintf('%s\n',filename.name)
    % Why? This works on a Mac?
    if ~strcmp(computer,'GLNXA64')
        filename=filename.name;
    end
    DATA{tr} = import_tsv_from_qtm_to_matlab(filename);
end
bodies_labels = {'violin1','violin2','viola','cello'};


%% Should we rotate the bodies so that the AP and ML axes correspond to Y 
% and X axes in the data frame, respectively?
target_vector = [-1 0];
dims = [1 2];
reference_markers = {'backl','backr'};
for tr = 1:numel(filenames)
    fprintf('%s\n',filenames(tr).name)
    DATA{tr} = rotate_to_target_vector(DATA{tr},target_vector,dims,bodies_labels,reference_markers,plotting_flag);
    if plotting_flag==1;pause;end
end

%% Get proportion of NaNs for each participant and marker.
% Trials are cells, markers are columns.
for tr = 1:numel(DATA)
    for marker=1:size(DATA{tr}.X,3)
        %DATA{tr}.missing_count(1,marker) = sum(isnan(DATA{tr}.X(:,1,marker))); % this is absolute count
        DATA{tr}.missing_prop(1,marker) = sum(isnan(DATA{tr}.X(:,1,marker))/size(DATA{tr}.X(:,1,marker),1)); %proportion
    end
end
% What to do with nan gaps? Decide depending on how broken is the data, 
% fill the nan parts or cut them out, or the whole marker?


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
                        fprintf('%6.3f,',DATA{tr}.missing_prop(m))
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
% We still have to verify that's the case for piece2.


%% To confirm the data looks OK set the flag to 1 and run an animated plot.
plotting_video_flag = 0;
if plotting_video_flag == 1
    for tr = 1:numel(DATA)
        fprintf('%s\n',filenames(tr).name)
        plot_animated_in_3d(DATA{tr}.X)
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

% Now prepare data for gc analysis - EW
% specify the names of the markers I want to analyze!
headMark = {'cellohat0','cellohat1','cellohat2','cellohat3','violahat0','violahat1','violahat2','violahat3',...
    'violin1hat0','violin1hat1','violin1hat2','violin1hat3','violin2hat0','violin2hat1','violin2hat2','violin2hat3'};

% Take the desired markers, downsample the data, convert to z-scores, and 
% take the average of the 4 head markers for each musician.
plotting_flag=0; %set this to 1 if you want to see the plots
DATA=prepare_data_for_mvgc(DATA,headMark,plotting_flag);

% Our DATA variable has a struct called AP now, which is the
% anterior-posterior body sway of the 4 musicians.

% Get the data into a matrix form for the MVGC toolbox
M=create_matrix_for_mvgc(DATA, plotting_flag);


%% What movement variables with reduced dimension should to use?
% v AP sway (front-back). (Y-axis in the rotated data).
% v MSD or speed. (Do that!)
% x PCA? (Probably not necessary here to fine-tune beyond ML and AP. We know that PC1 and PC2 will correspond closely to the AP and ML axes.)
% x The thing from quaternions? (No time for that.)
