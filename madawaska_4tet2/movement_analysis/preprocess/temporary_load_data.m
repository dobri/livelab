%% Where's my data?
if strcmp(computer,'GLNXA64')
    % '/media/dobri/disk2/mcmc_large_data/madawaska_tsv_data'
    data_folder = input('Where is my fokin money? '); % https://www.youtube.com/watch?v=pslUqRMPthk
else
    data_folder = '/Users/emilywood/Desktop/MATLAB/trainorlab/madawaska_4tet/data/ensemble';
end
which_piece = ['piece' num2str(input('Which of the two pieces? Type a number: '),'%1.0f')];
plotting_flag = input('To visualize data along the way? Type 1 for yes: ');

load(fullfile(data_folder,['DATA_' which_piece '.mat']));


% ...

dt=1/DATA{tr}.sf;a = [0;0;dot(diff(diff(X)./dt)./dt,diff(diff(X)./dt)./dt,2).^.5];a([1:300 end-299:end]) = 0;