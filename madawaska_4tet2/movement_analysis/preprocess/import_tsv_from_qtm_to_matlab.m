function out = import_tsv_from_qtm_to_matlab(filename)
%import_tsv_from_qtm_to_matlab(filename)
% This is supposed to be generic, but in practice there are always some
% subtle differences in how QTM exports its data files to a .tsv, so this
% script has to be edited for the given project.

if isempty(filename)
    'What am I supposed to do without input?';
    out = [];
    return
end


delimiter = '\t';
fileID = fopen(fullfile(filename.folder,filename.name),'r');
paramArray = textscan(fileID, '%s%f', 6, 'Delimiter', delimiter, ...
    'EmptyValue' ,NaN,'HeaderLines', 0, 'ReturnOnError', false, 'EndOfLine', '\r\n');


% Readoff some parameters.
for n=1:6
    if strcmp(paramArray{1,1}{n},'NO_OF_FRAMES')
        nrows=paramArray{1,2}(n);
    end
    if strcmp(paramArray{1,1}{n},'NO_OF_MARKERS')
        nchannels=paramArray{1,2}(n);
    end
    if strcmp(paramArray{1,1}{n},'FREQUENCY')
        sf=paramArray{1,2}(n);
    end
end


% Read the marker labels.
col_names = textscan(fileID, '%s', nchannels+1, 'Delimiter', delimiter, ...
    'EmptyValue' ,NaN,'HeaderLines', 4, 'ReturnOnError', false, 'EndOfLine', '\r\n');
col_names = col_names{1}(2:end);


% Read the data.
dataArray = textscan(fileID, '%f', 'Delimiter', delimiter, 'EmptyValue' ,NaN,...
    'HeaderLines', 2, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);


% Pull the raw numeric data. Then reshape it as an array with shape =
% number of samples x three dimensions x number of markers.
X=reshape(dataArray{1},nchannels*3,nrows)';
X=reshape(X,[],3,nchannels);
clear dataArray


out.filename = filename.name;
out.col_names = col_names';
out.X = X;
out.sf = sf;