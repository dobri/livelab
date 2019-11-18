function [X_detrended,The_trend] = preprocess(DATA,varargin)
if isempty(varargin{1})
    plotting_flag = 1;
else
    plotting_flag = varargin{1};
end

X = DATA.X;
sf = DATA.sf;
t = (1:size(X,1))'./sf; % create time variable


%% Plot each dimension against time for each marker.
if plotting_flag == 2
    for marker=1:size(X,3)
        for d=1:3
            plot(t(1:10:end),X(1:10:end,d,marker),'-')
            hold on
        end
        hold off
        xlabel('time (s)')
        ylabel('X,Y, or Z [mm]')
        legend('X','Y','Z')
        title([DATA.filename(~(double(DATA.filename)==95)) ' - ' DATA.col_names{marker}]) % ,'interpreter','latex'
        pause
    end
end


%% De-trend, de-nan, de-artefact.
X_detrended = zeros(size(X,1),size(X,2),size(X,3));
The_trend = nan(size(X,1),size(X,2),size(X,3));
for marker=1:size(X,3)
    [X_detrended(:,:,marker), The_trend(:,:,marker)] = remove_nonstation_and_nans_with_splines(X(:,:,marker),sf,t,[],[],plotting_flag);
    if plotting_flag == 1
        text(.45,-.3,[DATA.filename(~(double(DATA.filename)==95)) ' - ' DATA.col_names{marker}],'units','normalized')
        pause
    end
end


%% Get proportion of NaNs for each participant and marker.
% What to do with nan gaps? Decide depending on how broken is the data, 
% fill the nan parts or cut them out, or the whole marker?


%% Here, fill gaps by linearly interpolating, good enough for small gaps.
X_detrended([1 end],:,:)=0;
The_trend([1 end],:,:)=0;
for marker=1:size(X,3)
    for d=1:3
        for k=2:size(X,1)-1
            if isnan(X_detrended(k,d,marker))
                k_end = min([size(X_detrended,1) k-1+find(~isnan(X_detrended(k:end,d,marker)),1,'first')]);
                X_detrended(k:k_end-1,d,marker) = linspace(X_detrended(k-1,d,marker),X_detrended(k_end,d,marker),k_end-k);
                The_trend(k:k_end-1,d,marker) = linspace(The_trend(k-1,d,marker),The_trend(k_end,d,marker),k_end-k);
            end
        end
    end
end

