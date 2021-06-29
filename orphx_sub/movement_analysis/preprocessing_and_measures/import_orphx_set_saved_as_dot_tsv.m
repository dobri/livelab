function [D,sf,IDS_num,IDS_str] = import_orphx_set_saved_as_dot_tsv(filename,filename_markersids,plotting_flag)

delimiter = '\t';
fileID = fopen(filename,'r');
paramArray = textscan(fileID, '%s%f', 6, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', 0, 'ReturnOnError', false, 'EndOfLine', '\r\n');
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
col_names = textscan(fileID, '%s', nchannels+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', 4, 'ReturnOnError', false, 'EndOfLine', '\r\n');
col_names = str2double(col_names{1}(2:end));

% Pull the raw numeric data. Then reshape it as an array with shape =
% number of sample x three dimensions x number of markers.
dataArray = textscan(fileID, '%f', 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'HeaderLines', 1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
D=reshape(dataArray{1},nchannels*3,nrows)';
D=reshape(D,[],3,nchannels);
clear dataArray

% Check that we're getting the right number of individual markers.
if size(D,3)~=numel(col_names)
    fprint('F***!\n')
    keyboard
end

% Visually inspect the whole audience.
if plotting_flag == 1
    figure(1)
    if 0
        for pp=1:size(D,3)
            plot3(D(1:1e2:end,1,pp),D(1:1e2:end,2,pp),D(1:1e2:end,3,pp));
            text(nanmedian(D(1:1e2:end,1,pp)),nanmedian(D(1:1e2:end,2,pp)),nanmedian(D(1:1e2:end,3,pp))+50,num2str(col_names(pp),'%.0f'))
            hold on
        end
        grid on
        hold off
        xlabel('x')
        ylabel('y')
        zlabel('z')
    else
        % But it's easier in a 2D view from above.
        for pp=1:size(D,3)
            plot(D(1:1e2:end,1,pp),D(1:1e2:end,2,pp));
            text(nanmedian(D(1:1e2:end,1,pp)),nanmedian(D(1:1e2:end,2,pp)),num2str(col_names(pp),'%.0f'))
            hold on
        end
        grid on
        hold off
        xlabel('x')
        ylabel('y')
    end
end

% Labels. What do we mean by that? These were alphabet chars that allowed
% us to identify participants and link them anonymously to other modalities
% of measurement such as questionnaires and ratings.
% This will merge redundant markers and then return the first of the two on the head.
% Makers not associated with an alphanumeric label are not retained. I
% think there are a couple of these, presumably bad quality recording or left early.
[IDS_str,IDS_num] = import_subj_marker_labels(filename_markersids);
[D,IDS_str,IDS_num] = crop_D_and_labels(D,col_names,IDS_str,IDS_num); 

% Compare participant labels and marker numbers with those in Subj_Markers.txt
for r=1:numel(IDS_str)
    fprintf('%4s,%4.0f\n',IDS_str{r},IDS_num(r))
end

if plotting_flag == 1
    figure(2)
    for pp=1:size(D,3)
        plot(D(1:1e2:end,1,pp),D(1:1e2:end,2,pp));
        text(nanmedian(D(1:1e2:end,1,pp)),nanmedian(D(1:1e2:end,2,pp)),IDS_str{pp})
        hold on
    end
    grid on
    hold off
    xlabel('x')
    ylabel('y')
end

end


%%
function [D_merged,IDS_str_merged,IDS_num_merged] = crop_D_and_labels(D,col_names,IDS_str,IDS_num)

IDS_num_merged = zeros(size(IDS_str));

D_merged = zeros(size(D,1),3,size(IDS_num,2));
for n = 1:size(IDS_num,2)
    [r,~] = find(col_names==IDS_num{n});
    if 0 % debugging
        for d=1:3
            for rr=1:numel(r)
                plot(D(:,d,r(rr)))
                hold on
            end
            hold off
            pause
        end
    end
    % Here's the important tricky part. We merge by simply nan-summing,
    % which assumes that different vectors that stand for the same marker
    % will complement w/out overlapping each other, the empty parts are nans.
    % D_merged(:,:,n) = nansum(D(:,:,r),3);
    D_merged(:,:,n) = nanmean(D(:,:,r),3);
    for d=1:3
        D_merged(D_merged(:,d,n)==0,:,n) = nan;
    end
    IDS_num_merged(n) = col_names(r(1));
    if 0 % debugging
        disp(r)
        disp(IDS_num{n})
        disp(IDS_num_merged(n))
        disp(IDS_str{n})
        disp(col_names(r))
        
        plot(D_merged(:,:,n))
        hold on
        for rn=1:numel(r)
            plot(D(:,:,r(rn)),'--','linewidth',2)
        end
        hold off
        pause
    end
end

% Find the redundant vectors.
remove_markers = [];
for n = 2:size(IDS_num,2)
    % Assume that labels and matrices are ordered, so only seek contiguous repetitions.
    if strcmp(IDS_str{n-1},IDS_str{n})
        % Remove the marker with more nans.
        if sum(isnan(D_merged(:,1,n-1))) > sum(isnan(D_merged(:,1,n)))
            remove_markers = horzcat(remove_markers, n-1);
        else
            remove_markers = horzcat(remove_markers, n);
        end
    end
end

keep_markers = find(all((((1:numel(IDS_str))'==remove_markers)==0)'));
D_merged = D_merged(:,:,keep_markers);
IDS_num_merged = IDS_num_merged(keep_markers);
IDS_num_merged = floor(IDS_num_merged./10);
IDS_str_merged = IDS_str(keep_markers);

end


%%
function [IDS_str,IDS_num] = import_subj_marker_labels(filename)

delimiter = ',';
fileID = fopen(filename, 'r');
x = textscan(fileID, '%s', 'Delimiter', delimiter, 'EmptyValue', NaN, 'HeaderLines', 0, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

IDS_str = cell(1,1);
IDS_num = cell(1,1);
c = 0;
for n = 1:numel(x{1})
    if isnan(str2double(x{1}{n}))
        c = c+1;
        IDS_str{c} = x{1}{n};
        IDS_num{c} = [];
    else
        IDS_num{c} = [IDS_num{c} str2double(x{1}{n})];
    end
end
end