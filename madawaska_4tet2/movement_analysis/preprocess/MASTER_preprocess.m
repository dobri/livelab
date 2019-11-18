plotting_video_flag = 0;
plotting_flag = 0;
data_folder = 'D:\str4tet_madawaska\Data';
%cd(data_folder)
%filename = 'piece1_solo1.tsv';
filenames = dir(fullfile(data_folder,'*.tsv'));

for tr = 1%:numel(filenames)
    DATA{tr} = import_tsv_from_qtm_to_matlab(filenames(tr),plotting_video_flag);
end

for tr = 1:numel(DATA)
    [DATA{tr}.X_detrended,DATA{tr}.The_trend] = preprocess(DATA{tr},plotting_flag);
end