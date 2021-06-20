function [D_detrended_2,The_trend] = removal_of_nonstation_splines(D,sf,t,ds,smooth_param,plotting_flag)
%%removal_of_nonstation_splines
% D is nxdxk array, where n is number of samples, d is dimension (1:3), and
% k is number of markers, or participants, or whatever else is the unit.

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

    t_index = find(sum(isnan(x),2)==0);
    t_prim = (1:numel(t_index))';
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
            subplot(4,1,1)
            plot(t,D(:,d,pp))
            
            subplot(4,1,2)
            plot(t_prim./sf,x(:,d),t_prim./sf,the_trend(:,d))
            
            subplot(4,1,3)
            plot(t_prim./sf,fitted)
            
            subplot(4,1,4)
            plot(t,x_stationary(:,d))
            pause
        end
    end
end