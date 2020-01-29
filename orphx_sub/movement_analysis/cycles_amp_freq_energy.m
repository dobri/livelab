function cycl = cycles_amp_freq_energy(x,t,sf,varargin)
% Remove outliers.
% Interpolate the amplitude time-series.
% For H energy, set the x where the amps < 20 mm as 0.
% But not the v. In this way, the H will be at least v^2, plus
% potential energy if any.

if isempty(varargin)
    plotting = 1;
else
    plotting = varargin{1};
end


%% Peak-picking.
min_peak_prom = .2;
min_period = ceil(.333*sf);
for d=1:size(x,2)
    [pks1{d},locs1{d}]=findpeaks( x(:,d),'MinPeakDistance',min_period,'MinPeakProminence',min_peak_prom);
    [pks2{d},locs2{d}]=findpeaks(-x(:,d),'MinPeakDistance',min_period,'MinPeakProminence',min_peak_prom);
end


%% Break traj in cycles and get amplitudes and periods at the peaks.
cycles=cell(1,size(x,2));
for d=1:size(x,2)
    cycles{d}=zeros(size(pks1{d},1)-1,size(x,2));
    for cc=1:size(cycles{d},1)
        neg_peak = find((t(locs2{d}) > t(locs1{d}(cc))).*(t(locs2{d}) < t(locs1{d}(cc+1))), 1);
        if ~isempty(neg_peak)
            cycles{d}(cc,1)=range(x(locs1{d}(cc):locs1{d}(cc+1),d))/2; % half-amplitudes
            cycles{d}(cc,2)=t(locs1{d}(cc+1))-t(locs1{d}(cc)); % periods
        end
    end
    cycles{d}(:,3)=cycles{d}(:,1).^2.*(2*pi./cycles{d}(:,2)).^2./2; % energies
end


%% Clean outlier cycles.
for d=1:size(x,2)
    index_remove{d}=false(size(cycles{d},1),1);
    for dv=1:size(cycles,2)
        index_remove{d}(cycles{d}(:,dv) > (median(cycles{d}(:,dv))+4*std(cycles{d}(:,dv))))=true;
        index_remove{d}(isnan(cycles{d}(:,dv)))=true;
    end
    index_remove{d}(cycles{d}(:,1)<10)=true;
end


%% Interpolate the period from peak to peak to make a time-series.
period_as_ts=zeros(size(x));
amplit_as_ts=zeros(size(x));
for d=1:size(x,2)
    if ~isempty(cycles{d})
        period_as_ts(1:locs1{d}(2),d)=nanmean(cycles{d}(1:2,2));
        amplit_as_ts(1:locs1{d}(2),d)=nanmean(cycles{d}(1:2,1));
        for cc=3:size(cycles{d},1)
            if index_remove{d}(cc-1)~=1
                %                 period_as_ts(locs1{d}(cc-1):locs1{d}(cc),d)=cycles{d}(cc-2,2);
                %             else
                %                 period_as_ts(locs1{d}(cc-1):locs1{d}(cc),d)=cycles{d}(cc-1,2);
                %             end
                period_as_ts(locs1{d}(cc-1):locs1{d}(cc),d)=cycles{d}(cc-1,2);
                amplit_as_ts(locs1{d}(cc-1):locs1{d}(cc),d)=cycles{d}(cc-1,1);
            end
        end
        period_as_ts(locs1{d}(end):end,d)=cycles{d}(end,2);
        amplit_as_ts(locs1{d}(end):end,d)=cycles{d}(end,1);
    end
    %keyboard
    %period_as_ts(period_as_ts(:,d)==0,d)=median(period_as_ts(:,d));
    cycles{d}(index_remove{d},:)=[];
end


% Just add the durations of all non index_remove rows in cycle. Then divide
% by the total trial lenght. No, because these periods in PC1 and PC2 may 
% or may not overlap.
oscillatory = zeros(size(x));
for d=1:size(x,2)
    for cc=2:size(cycles{d},1)
        if index_remove{d}(cc)~=1
            oscillatory(locs1{d}(cc-1):locs1{d}(cc),d)=1;
        end
    end
end
oscillatory = sum(oscillatory,2)>0;
osc_prop = sum(oscillatory)/size(t,1);


%% H
m=80;
H=zeros(size(x,1),size(x,2)+1);
Q=zeros(size(x,1),size(x,2));
q=zeros(size(x,1),2);
for d=1:size(x,2) % No point of doing the ML and AP sway. Not sure what it's related to.
    vel=[0;diff(x(:,d))./(1/sf)];
    q(:,2)=m*vel.^2;
    q(:,1)=(2*pi./period_as_ts(:,d)).^2.*(x(:,d)).^2;
    q(q(:,1)==inf,1)=0;
    q(q(:,2)==inf,2)=0;
    H(:,d)=((q(:,2)+q(:,1))./2).^.5;
    Q=Q+q;
end
H(:,size(x,2)+1)=sum(H(:,1:size(x,2)),2);


%% Diagnostic plot.
if plotting>0
    figure(2)
    subplot(3,1,1)
    plot(t,x,'linewidth',2);
    for d=1:size(x,2)
        hold on
        plot(t(locs1{d}),x(locs1{d},d),'v','linewidth',2)
        plot(t(locs2{d}),x(locs2{d},d),'^','linewidth',2)
    end
    hold off

    subplot(3,1,2)
    plot(t,period_as_ts,'linewidth',2)
    hold on
    plot(t,amplit_as_ts,'linewidth',2)
    plot(t,oscillatory*10,'linewidth',2)
    hold off

    subplot(3,1,3)
    plot(t,H,'linewidth',2)
    
    if plotting == 1
        pause
    end
end


cycl.osc_prop = osc_prop;
cycl.cycles = cycles; % iois in column 2.
cycl.peak_loc = locs1;
cycl.peak_amp = pks1;
cycl.amp_interp = amplit_as_ts;
cycl.freq_interp = 1./period_as_ts;
cycl.h = H;
cycl.q = Q;
