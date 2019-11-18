function [D_detrended_2,The_trend] = remove_nonstation_and_nans_with_splines(D,sf,t,ds,smooth_param,plotting_flag)
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


D_detrended_2 = zeros(size(D,1),size(D,2),size(D,3));
The_trend = nan(size(D,1),size(D,2),size(D,3));
for pp=1:size(D,3)
    
    x = D(:,:,pp);
    
    % prep time and x vectors w/out nan values.
    t_index = find(sum(isnan(x),2)==0);
    t_prim = (1:numel(t_index))';
    t_prim_plotting = t(sum(isnan(x),2)==0);
    x = x(sum(isnan(x),2)==0,:);
    x_stationary = nan(size(D,1),3);
    
    the_trend = zeros(size(x));
    for d=1:3
        
        f=fit(t_prim(1:ds:end),x(1:ds:end,d),'smoothingspline','SmoothingParam',smooth_param);
        the_trend(:,d) = feval(f,t_prim);
        
        fitted = x(:,d)-the_trend(:,d);
        
        for n=1:numel(fitted)
            x_stationary(t_index(n),d) = fitted(n);
            The_trend(t_index(n),:,pp) = the_trend(n,:);
        end
        D_detrended_2(:,:,pp) = x_stationary;
        
        if plotting_flag == 1
            subplot(6,1,(d-1)*2+1)
            plot(t,D(:,d,pp))
            hold on
            %plot(t_prim./sf,x(:,d),'--')
            plot(t_prim_plotting,x(:,d),'--')
            plot(t_prim_plotting,the_trend(:,d))
            hold off
            ylabel([char(87+d) ', mm'])
            if d==1
                legend('data','data w/out nans','the trend')
            end
            
            subplot(6,1,(d-1)*2+2)
            plot(t_prim_plotting,fitted)
            hold on
            plot(t,x_stationary(:,d),'--')
            plot(t,isnan(D(:,d,pp))*2e2)
            hold off
            ylabel([char(87+d) ', mm'])
            if d==1
                legend('detrended data, w/out nans','... w/ nans','nans')
            end
        end
    end
    
    %if plotting_flag == 1
    %    pause
    %end
end