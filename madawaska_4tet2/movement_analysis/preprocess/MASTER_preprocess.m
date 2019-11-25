%% Where's my data?
plotting_video_flag = 0;
plotting_flag = 0;
data_folder = input('Type the path to your .tsv data files, i.e. ''D:\\str4tet_madawaska\\Data\'' or ''/Users/emilywood/Desktop/MATLAB/madawaska_4tet/data''\n>> ');
fprintf('%s\n',['We''ll be looking inside ' data_folder])
filenames = dir(fullfile(data_folder,'*.tsv'));


%% Import the QTM data saved as text files.
for tr = 1:numel(filenames)
    filename = filenames(tr); %filename = 'piece1_solo1.tsv'
    fprintf('%s\n',filename.name)
    DATA{tr} = import_tsv_from_qtm_to_matlab(filename);
end


%% Should we rotate the bodies so that the AP and ML axes correspond to Y 
% and X axes in the data frame, respectively?
target_vector = [-1 0];
dims = [1 2];
reference_markers = {'backl','backr'};
bodies_labels = {'violin1','violin2','viola','cello'};
for tr = 1:numel(filenames)
    fprintf('%s\n',filenames(tr).name)
    DATA{tr} = rotate_to_target_vector(DATA{tr},target_vector,dims,bodies_labels,reference_markers);
    pause
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
wanted_head_markers{1} = {'violin1hat0','violin1hat1','violin1hat2','violin1hat3'};
wanted_head_markers{2} = {'violin2hat0','violin2hat1','violin2hat2','violin2hat3'};
wanted_head_markers{3} = {'violahat0','violahat1','violahat2','violahat3'};
wanted_head_markers{4} = {'cellohat0','cellohat1','cellohat2','cellohat3'};

for tr = 1:numel(DATA)
    for body = 1:4
        for marker = 1:numel(wanted_head_markers{body})
            for m = 1:numel(DATA{tr}.col_names)
                if strcmp(DATA{tr}.col_names{m},wanted_head_markers{body}{marker})
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
% It seems that head markers were perfectly recorded, at least for piece1.
% We still have to verify that's the case for piece2.


%% To confirm the data looks OK set the flag to 1 and run an animated plot.
if plotting_video_flag == 1
    for tr = 1:numel(DATA)
        fprintf('%s\n',filenames(tr).name)
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
    DATA{tr}.the_Trend = nan(size(DATA{tr}.X,1),size(DATA{tr}.X,2),size(DATA{tr}.X,3));
    for marker=1:size(DATA{tr}.X,3)
        [DATA{tr}.X_detrended(:,:,marker), DATA{tr}.the_Trend(:,:,marker)] = remove_nonstation_and_nans_with_splines(DATA{tr}.X(:,:,marker),DATA{tr}.sf,(1:size(DATA{tr}.X,1))'./DATA{tr}.sf,[],[],plotting_flag);
        if plotting_flag == 1
            text(.45,-.3,[DATA{tr}.filename(~(double(DATA{tr}.filename)==95)) ' - ' DATA{tr}.col_names{marker}],'units','normalized')
            pause
        end
    end
end


%% What movement variables with reduced dimension should we use?
% v AP sway (front-back). (Y-axis in the rotated data).
% x PCA? (Probably not necessary here to fine-tune beyond ML and AP.)
% -> MSD or velocity (Do that!)
% x The thing from quaternions? (No time for that.)
