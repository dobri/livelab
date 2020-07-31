%% Preprocess data for MVGC toolbox
% Lucas Klein
% June 2020

% This script reads in data from each participant stored as text files and
% assigns that data to variables with the following convention: 



%% Txt files input

clear all
clc
cd '~/Desktop/Musical_following/ANALYSIS/Txts'
path = '~/Desktop/Musical_following/ANALYSIS/Txts';
path_stim = '~/Desktop/Musical_following/ANALYSIS/Txts_stim';
data_folder = dir();



%% PARAMETERS

sf = 44100;

num_of_participants = 1;

ds_targets = [10, 50, 100, 200, 300, 400, 500, 1000, 1500, 2000];

plotting_flag = 1;



%% Loop through folders in path to make a list of participants
participants = string(zeros(1,num_of_participants)); % empty list

for p = 4:length(data_folder) % Because the 1st, 2nd and 3rd elements in data_folder.name are nonsense
    par = data_folder(p).name;
    participants(p-3) = par; % fill the participant list
end
% participants is now an arrray of strings with all participants' names



%% Loop through participants
% First, loop through each participant in data_folder.
for p = 1:numel(participants)
    participant_folder = append(path,'/',participants(p));
    participant_folder_stim = append(path_stim,'/',participants(p));
    addpath(participant_folder);
    addpath(participant_folder_stim);
    filenames = dir(fullfile(participant_folder,'*.txt'));
    filenames_stim = dir(fullfile(participant_folder_stim,'*.txt'));
    participant_data = cell(1,numel(filenames));
    participant_data_stim = cell(1,numel(filenames_stim));
    
    
    %% Load the stimulus track
%     stim = load('~/Desktop/Musical_following/ANALYSIS/arp1.txt'); % Load the recording
%     stim = stim(:,1); % Which column to use? (1 or 2)


    
    %% Now load the performances
    % Loop through all files in each participant_folder (trials)
    % Within consecutive cells in DATA (cell array), make a new struct for each
    % trial and load data into a field within each struct called X
    
    
    
    for trial = 1:numel(filenames) % change and repeat for stim
        filename = filenames(trial).name;
        filename_stim = filenames_stim(trial).name;
        X = load(filename);
        stim = load(filename_stim);
        X = X(:,1); % Which column to use? (1 or 2)
        stim = X(:,1); % Which column to use? (1 or 2)
        
        % Truncate both time series to the length of the shortest
        
        trunc = min(length(X),length(stim));
        X = X(1:trunc,1);
        stim = stim(1:trunc,1);
        
        %% Loop through all the downsampling rates we want to check
        % For each one, downsample and z-score both the stimulus and
        % performance using the prepare_following function
        
        for ds_target = 100
            
            % Downsample and z-score the data
            stim_prep = prepare_following(stim, ds_target, sf, plotting_flag);
            X_prep = prepare_following(X, ds_target, sf, plotting_flag);

            combined = cat(3,X_prep,stim_prep);
            
            ds_label = ['ds_prep_' + string(ds_target)];
            DATA{trial}.(ds_label) = combined; % Make new field for each trial
            DATA{trial}.participant = participants(p);
            %DATA{trial}.fields = fieldnames(DATA{trial});
            
            % At this point, DATA has cells for each trial, each of which contain
            % one struct for each downsample target (and one for participant name)            
            
            % Optional plotting
            
            
            %% SAVE THE DATA
            % Save a NEW variable D that contains a cell for each participant, each of
            % which contains one matrix for each downsample target

            plotting_flag = 0; % Set to 1 if you want to plot all the matrices

            M_label = ['M_' + string(ds_target)];
            D{p}.(M_label) = create_matrix_following(DATA, plotting_flag);

        end


    end
    
    
end


% Now take this variable D over to the mvgc toolbox for gc!
save('/Users/lucas/Desktop/Musical_following/ANALYSIS/D','D')



