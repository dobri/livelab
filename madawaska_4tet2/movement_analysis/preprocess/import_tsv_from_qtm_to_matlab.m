function out = import_tsv_from_qtm_to_matlab(varargin)
if isempty(varargin{1})
    filename = 'piece1_solo1.tsv';
else
    filename = varargin{1};
end
if isempty(varargin{2})
    plotting_video_flag = 1;
else
    plotting_video_flag = varargin{2};
end
%cd(data_folder)


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


%% Video-like plotting.
if plotting_video_flag == 1
    x1=squeeze(X(1,1,:));
    y1=squeeze(X(1,2,:));
    z1=squeeze(X(1,3,:));
    p = plot3(x1,y1,z1,'ko','MarkerSize',5,'MarkerFaceColor','r');
    %p.XDataSource = 'x1';
    %p.YDataSource = 'y1';
    %p.ZDataSource = 'z1';
    grid on
    xlim([-1000 3000])
    ylim([-0 3200])
    zlim([  500 1500])
    set(gca,'view',[-13 74])
    for t=2:size(X,1)
        x1=squeeze(X(t,1,:));
        y1=squeeze(X(t,2,:));
        z1=squeeze(X(t,3,:));
        %refreshdata
        set(p,'XData',x1)
        set(p,'YData',y1)
        set(p,'ZData',z1)
        drawnow
        pause(1/sf)%pause
    end
end

out.filename = filename.name;
out.col_names = col_names;
out.X = X;
out.sf = sf;