%% Use this script (not MASTER_preprocess) IF YOU ALREADY HAVE DATA FILES FROM DOBRI (DATA_piece1.mat, DATA_piece2.mat)

% first, we'll define some variables 
bodies_labels = {'violin1','violin2','viola','cello'};
headMark = {'cellohat0','cellohat1','cellohat2','cellohat3','violahat0','violahat1','violahat2','violahat3',...
    'violin1hat0','violin1hat1','violin1hat2','violin1hat3','violin2hat0','violin2hat1','violin2hat2','violin2hat3'};
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
    %DATA2{tr}.X_clean=DATA2{tr}.X;
    DATA2{tr}.X_clean = zeros(size(DATA2{tr}.X,1),size(DATA2{tr}.X,2),size(DATA2{tr}.X,3));
     for marker=1:size(DATA2{tr}.X,3)
        DATA2{tr}.X_clean(:,:,marker) = fill_nans_with_splines_X(DATA2{tr}.X(:,:,marker));
    end
end

%check to see if it worked
for tr = 1:numel(DATA2)
    fprintf('Trial %1.0f\n',tr)
    for marker=1:size(DATA2{tr}.X_clean,3)
        if any(strcmp(headMark,DATA2{tr}.col_names{marker}))
            fprintf('%6.4f\n',sum(isnan(DATA2{tr}.X_clean(:,1,marker)))/size(DATA2{tr}.X_clean(:,1,marker),1)); %proportion of nans. should all be 0
        end
    end
end

%Add an X_clean into DATA1 to keep things parallel
for i = 1:numel(DATA1)
    DATA1{i}.X_clean=DATA1{i}.X;
end

%DATA1=DATA(1:8);
%DATA2=DATA(9:16);

%% Now prepare data for gc analysis 

% Take our two data trajectories 'X' and 'X_detrended' and then take 
% our desired markers (head markers), downsample the data, convert to  
% z-scores, and take the average of the 4 head markers for each musician.

dataTrajs={'X_clean','X_detrended'};
for i = 1:length(dataTrajs)
    DATA1=prepare_data_for_mvgc(DATA1,dataTrajs{i},headMark,plotting_flag);
    DATA2=prepare_data_for_mvgc(DATA2,dataTrajs{i},headMark,plotting_flag);
end

% Our DATA variables have new structs called X_processed and
% X_detrended_processed, which correspond to the anterior-posterior body 
% sway of the 4 musicians in each of these trajectories.

% Get the data into a matrix form for the MVGC toolbox
dataTrajs={'X_clean_processed','X_detrended_processed'};
for i = 1:length(dataTrajs)
  D{1}.(dataTrajs{i})=create_matrix_for_mvgc(DATA1, dataTrajs{i}, plotting_flag);
  D{2}.(dataTrajs{i})=create_matrix_for_mvgc(DATA2, dataTrajs{i}, plotting_flag);
end

% Now we have a NEW data variable called D.
% There are TWO structs in it corresponding to piece 1 and piece 2
% WIthin these structs, we have a matrix for mvgc for each of our two
% data trajectories, including X and X_detrended


% Add acceleration data into 'D'
for i = 1:length(DATA1)
    A1{i}=permute(DATA1{i}.A, [2,1]);
    D{1}.A=A1;
    A2{i}=permute(DATA2{i}.A, [2,1]);
    D{2}.A=A2;
end

%Check acceleration data. I noticed there are a lot of 0s at the end of most trials...
%plot(D{1}.A{1}(1,:))
%This is fine

%Check data - make sure there is no trend that will affect CC analysis
if plotting_flag ==1
    for piecei = 1:2
        for triali=1:8 
            figure
            cel=D{piecei}.X_clean_processed(1,:,triali); 
            vl=D{piecei}.X_clean_processed(2,:,triali);
            v1=D{piecei}.X_clean_processed(3,:,triali);
            v2=D{piecei}.X_clean_processed(4,:,triali);
            plot(cel)
            hold on
            plot(vl)
            hold on
            plot(v1)
            hold on
            plot(v2)
            pause
        end
    end
end


% Now take this variable D over to the mvgc toolbox for gc!
%save('D','D')
