clear

% Move to where the data is stored, start reading-in parameters.
% i.e. '/media/dobri/disk2/mcmc_large_data/orphx'
data_path = input('Me: Where is my money?\nYou: ');
cd(data_path)
filename = fullfile(data_path,'orphx_set0005.tsv');
filename_markers_ids = fullfile(data_path,'Subj_Markers.txt');

% Set this to visually inspect the multiple stages of data processing.
plotting_flag = 1;

[D,sf,IDS_num,IDS_str] = import_orphx_set_saved_as_dot_tsv(filename,filename_markers_ids,plotting_flag);
t = (1:size(D,1))'./sf; % create time variable


%% Plot each dimension against time for each marker.
if plotting_flag == 1
    close all
    for pp=1:size(D,3)
        for d=1:3
            plot(t,D(:,d,pp))
            hold on
        end
        hold off
        xlabel('time (s)')
        ylabel('X,Y, or Z [mm]')
        legend('X','Y','Z')
        title(['Marker ',num2str(IDS_num(pp)) ', pp ' IDS_str{pp}])
        pause
    end
end


%% Plot 2D+time movement of each marker. (Not really useful.)
if 0
    for pp=1:size(D,3)
        plot3(t,D(:,1,pp),D(:,2,pp))
        xlabel('time (s)')
        ylabel('X (mm)')
        zlabel('Y (mm)')
        title(['Marker ',num2str(IDS_num(pp))])
        grid on
        ylim([-6e3 4e3])
        zlim([-4e3 1e4])
        pause
    end
end


%% Compute RMS for each marker for each dimension and place in an array.
% Subtract the mean from each score before computing RMS.
rms_array = zeros(size(D,3), 4);
for pp=1:size(D,3)
    for d=1:3
        rms_array(pp,d) = rms((D(:,d,pp))-mean(D(:,d,pp),'omitnan'),'omitnan');
    end
    %Add one more column to the rms array for the full-dimension rms.
    rms_array(pp,4) = rms((dot(D(:,:,pp)-nanmean(D(:,:,pp)),D(:,:,pp)-nanmean(D(:,:,pp)),2).^.5),'omitnan');
end


%% Windowed detrending using two different methods. Visually inspect.
close all
% [D_detrended_1, The_trend_1] = removal_of_nonstation_piecewise_lin_detrend(D,sf,t);
[D_detrended_2, The_trend_2] = removal_of_nonstation_splines(D,sf,t,[],[],plotting_flag);

% This variable manually codes the artifacts remaining after detrend and clean operations.
% Manually mark the points, such as when participants change spatial location by moving to a different spot in the room.
artifacts = {[],[],[],[],[],[],[207,213],[],[],[],[],[],[2923],[2725,2933],[],[733,815,3284,3285],[],[],[2506,2808],[],[],[],[],[2019,2154,2523],[],[1751,2677],[323,2115],[],[1088,1182,2109,2126],[19,1078],[1212,1952,2329],[],[],[270,1070,1370,1437,1706,2114,2195,2209,2959],[941],[1654],[272],[],[],[],[],[],[1635,2032,2599]};

% clear 2 seconds before and after a marked artifact point.
for pp=1:size(D,3)
    for d=1:3
        for a=1:numel(artifacts{pp})
            D_detrended_2(((artifacts{pp}(a)-2)*sf):((artifacts{pp}(a)+2)*sf),d,pp)=nan;
        end
    end
end


if plotting_flag == 1
    close all
    for pp=1:size(D,3)
        for d=1:3
            subplot(2,1,1)
            plot(t,D(:,d,pp))
            hold on
            % plot(t,The_trend_1(:,d,pp))
            plot(t,The_trend_2(:,d,pp))
            hold off
            
            subplot(2,1,2)
            plot(t,D_detrended_2(:,d,pp))
            
            hold on
            for a=1:numel(artifacts{pp})
                plot(artifacts{pp}(a),nanmean(D_detrended_2(:,d,pp)),'^','markersize',10,'linewidth',3)
            end
            plot(3350,nanmean(D_detrended_2(:,d,pp)),'<','markersize',10,'linewidth',3)
            hold off
            
            fprintf('%3.0f,%3.0f\n',pp,d)
            pause
        end
    end
end


%% Section the detrended data based on trial
trial_cut_points = [0.01 304;305 453;454 607;608 751;752 903;904 1051;1052 1199;1200 1349;1350 1499;1500 1649;1650 1800;1801 1951;1952 2099;2100 2249;2250 2399;2400 2550;2551 2701;2702 3350];
num_trial = size(trial_cut_points,1);
index = trial_cut_points*sf; % Trial beginning and end points in seconds, index = time * sf

for j = 1:length(index)
    D_trials{j} = D_detrended_2(index(j,1):index(j,2),:,:);
end


%% Get proportion of NaNs for each participant
PropNaN = zeros(num_trial,size(D_trials{1},3)); %43 participants x 18 trials
for pp = 1:length(PropNaN)
    for tr = 1:size(PropNaN,1)
        PropNaN(tr,pp)= sum(isnan(D_trials{tr}(:,1,pp)))/numel(D_trials{tr}(:,1,pp));
    end
end

% What should be the exclusion rule?
% Number of trials with at least half of the data missing.
nmissing(:,1) = sum(PropNaN >0.5)./num_trial;
% Number of full trials missing.
nmissing(:,2) = sum(PropNaN==1.0)./num_trial;


%% This is what will be used in later stages of the analysis.
ll.D = D_trials';
ll.IDSnum = IDS_num;
ll.IDSstr = IDS_str;
ll.nanprop = PropNaN;
ll.nan_trial_counts = nmissing';
ll.num_trial = num_trial;
ll.num_pp = size(D_trials{1},3);
ll.rms_raw_data = rms_array';
ll.sf = sf;
ll.trial_cut_points = trial_cut_points;


%{
 Save all of this into a mat file.
 save(fullfile(data_path,['orphx_set2_preprocessed_and_cleaned_',datestr(now,'yyyy-mm-dd'),'.mat']),'ll')
%}


%% That's it.