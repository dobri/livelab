%% Use this script (not MASTER_preprocess) IF YOU ALREADY HAVE DATA FILES FROM DOBRI (DATA_piece1.mat, DATA_piece2.mat)

% first, we'll define some variables 
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
%I already know that DATA1 markers were recorded perfectly. So let's check DATA2
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

% redo speed of head movement and calculate acceleration
DATA=[DATA1,DATA2];
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
        
        %Acceleration
        dt=1/DATA{tr}.sf;
        a = [0;0;dot(diff(diff(X)./dt)./dt,diff(diff(X)./dt)./dt,2).^.5];
        a([1:300 end-299:end]) = 0;
        DATA{tr}.A(:,b) = a;
         
        % Verify that there aren't spikes and nans remaining by accident.
        % Did we do a good job cleaning and filtering without killing v?
        if plotting_flag==1
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
end

DATA1=DATA(1:8);
DATA2=DATA(9:16);

%% Now prepare data for gc analysis 
% specify the names of the markers I want to analyze
headMark = {'cellohat0','cellohat1','cellohat2','cellohat3','violahat0','violahat1','violahat2','violahat3',...
    'violin1hat0','violin1hat1','violin1hat2','violin1hat3','violin2hat0','violin2hat1','violin2hat2','violin2hat3'};

% Take our two data trajectories 'X' and 'X_detrended' and then take 
% our desired markers (head markers), downsample the data, convert to  
% z-scores, and take the average of the 4 head markers for each musician.

dataTrajs={'X','X_detrended'};
for i = 1:length(dataTrajs)
    DATA1=prepare_data_for_mvgc(DATA1,dataTrajs{i},headMark,plotting_flag);
    DATA2=prepare_data_for_mvgc(DATA2,dataTrajs{i},headMark,plotting_flag);
end

% Our DATA variables have new structs called X_processed and
% X_detrended_processed, which correspond to the anterior-posterior body 
% sway of the 4 musicians in each of these trajectories.

% Get the data into a matrix form for the MVGC toolbox
dataTrajs={'X_processed','X_detrended_processed'};
for i = 1:length(dataTrajs)
  D{1}.(dataTrajs{i})=create_matrix_for_mvgc(DATA1, dataTrajs{i}, plotting_flag);
  D{2}.(dataTrajs{i})=create_matrix_for_mvgc(DATA2, dataTrajs{i}, plotting_flag);
end

% Now we have a NEW data variable called D.
% There are TWO structs in it corresponding to piece 1 and piece 2
% WIthin these structs, we have a matrix for mvgc for each of our two
% data trajectories, including X and X_detrended


% Now take this variable D over to the mvgc toolbox for gc!
save('D','D')
