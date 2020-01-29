function [D_detrended,The_trend] = remove_nonstation_and_nans_with_splines(D,sf,t,ds,smooth_param,plotting_flag,filename,col_names)
%%removal_of_nonstation_splines
% D is nxdxk array, where n is number of samples, d is dimension (1:3), and
% k is number of markers, or participants, or whatever else is the unit.
% start with addpath('/home/dobri/...')
% and
% [D,sf,IDS_num,IDS_str] = import_orphx_set_saved_as_dot_tsv();

if isempty(ds)
    ds=1e2;
end
if isempty(smooth_param)
    smooth_param = 1e-6;
end
if isempty(t)
    t = (1:size(D,1))'./sf;
end
if isempty(plotting_flag)
    plotting_flag = 1;
end


D_detrended = zeros(size(D,1),size(D,2),size(D,3));
The_trend = nan(size(D,1),size(D,2),size(D,3));
for pp=1:size(D,3)
    x = D(:,:,pp);
    
    % Prep time and x vectors w/out nan values.
    t_index = find(sum(isnan(x),2)==0);
    t_prim = (1:numel(t_index))';
    t_prim_plotting = t(sum(isnan(x),2)==0);
    x = x(sum(isnan(x),2)==0,:);
    x_stationary = nan(size(D,1),3);
    the_trend = zeros(size(D,1),3);
    for d=1:3
        
        f=fit(t_prim(1:ds:end),x(1:ds:end,d),'smoothingspline','SmoothingParam',smooth_param);
        trend_temp = feval(f,t_prim);
        fitted = x(:,d) - trend_temp;

        for n=1:numel(fitted)
            x_stationary(t_index(n),d) = fitted(n);
            the_trend(t_index(n),d) = trend_temp(n);
        end
        
        % Fill small nan gaps by linearly interpolating. Be careful with larger gaps!
        x_stationary([1 end],d)=0;
        the_trend([1 end],d)=nanmean(the_trend([1 end],d));
        for k=2:size(x_stationary,1)-1
            if isnan(x_stationary(k,d))
                k_end = min([size(x_stationary,1) k-1+find(~isnan(x_stationary(k:end,d)),1,'first')]);
                x_stationary(k:k_end-1,d) = linspace(x_stationary(k-1,d),x_stationary(k_end,d),k_end-k);
                the_trend(k:k_end-1,d) = linspace(the_trend(k-1,d),the_trend(k_end,d),k_end-k);
            end
        end
        
        % Visually inspect the data.
        if plotting_flag == 1
            subplot(6,1,(d-1)*2+1)
            plot(t,D(:,d,pp))
            hold on
            plot(t_prim_plotting,x(:,d),'--')
            plot(t,the_trend(:,d))
            hold off
            ylabel([char(87+d) ', mm'])
            if d==1
                legend('data','data w/out nans','the trend')
            end
            
            subplot(6,1,(d-1)*2+2)
            plot(t_prim_plotting,fitted)
            hold on
            plot(t,x_stationary(:,d),'--')
            % Here we show the nans.
            plot(t,isnan(D(:,d,pp))*2e2)
            hold off
            ylabel([char(87+d) ', mm'])
            if d==1
                legend('detrended data, w/ nans','... w/out nans','nans')
            end
        end
    end
    D_detrended(:,:,pp) = x_stationary;
    The_trend(:,:,pp) = the_trend;
    
    % Plot each dimension against time for each marker. Rarely useful.
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
            if ~isempty(filename) && ~isempty(col_names)
                title([filename(~(double(filename)==95)) ' - ' col_names{marker}]) % ,'interpreter','latex'
            end
            pause
        end
    end
end