%% Preprocess data for MVGC toolbox
% Lucas Klein
% June 2020

% This script reads in data from each participant stored as text files and
% assigns that data to variables with the following convention: 



%% Txt files input

clc
clear all
cd '~/Desktop/Musical_following/ANALYSIS/Txts'

path = '~/Desktop/Musical_following/ANALYSIS/Txts';
data_folder = dir();



%% Loop through folders in path to make a list of participants
num_of_participants = 3;
participants = string(zeros(1,num_of_participants)); % empty list

for p = 4 : length(data_folder) % Because the 1st, 2nd and 3rd elements in data_folder.name are nonsense
    par = data_folder(p).name;
    participants(p-3) = par; % fill the participant list
end
% participants is now an arrray of strings with all participants' names



%% Loop through participants
% First, loop through each participant in data_folder.
for p = 1 %:numel(participants)
    participant_folder = append(path,'/',participants(p));
    addpath(participant_folder);
    filenames = dir(fullfile(participant_folder,'*.txt'));
    participant_data = cell(1,numel(filenames));
    
    
    %% Load the stimulus track
    
    stim = load('~/Desktop/Musical_following/ANALYSIS/arp1.txt'); % Load the recording
    stim = stim(:,1); % Which column to use? (1 or 2)

    
    %% Now load the performances
    % Loop through all files in each participant_folder (trials)
    % Within consecutive cells in DATA (cell array), make a new struct for each
    % trial and load data into a field within each struct called X
    
    for trial = 1:numel(filenames)
        filename = filenames(trial).name;
        X = load(filename);
        X = X(:,1); % Which column to use? (1 or 2)
        
        % Truncate both time series to the length of the shortest
        
        trunc = min(length(X),length(stim));
        stim = stim(1:trunc,1);
        X = X(1:trunc,1);
        
        %% Downsample the data
        
        for ds_target = 100:100:2000
            stim = downsample_following(stim, ds_target);
            X = downsample_following(X, ds_target);
            
            combined = cat(3,X,stim);
            
            label = ['ds_' + string(ds_target)];
            DATA{trial}.(label) = combined; % Make new field for each trial
            DATA{trial}.participant = participants(p);
            %DATA{trial}.fields = fieldnames(DATA{trial});
            
        end
        
%         combined = cat(3,X,stim);
%         DATA{trial}.X = combined;
%         DATA{trial}.participant = participants(p);


    end
    
    
    
    %% Plot the data?
    plotting_flag = 0; % Set to 1 if you want to plot both matrices
    
    %% Save the data
    % Save a variable D that contains a cell for each participant's data
    % matrix M
    D{1}.M = create_matrix_following(DATA, plotting_flag);
    
    % D is a NEW matrix that contains one struct for each participant
    % Within each participant's struct, there is a matrix for each of the
    % downsampling rates we want to try (in progress)
    
end


% Now take this variable D over to the mvgc toolbox for gc!
save('/Users/lucas/Desktop/Musical_following/ANALYSIS/D','D')



