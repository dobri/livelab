function tempo_stats = tempo_modes_with_adaptive_kde(x,varargin)
% We need an adaptive method for kernel density estimation. So we adapt the
% bandwidth (in kde), the density of the grid (n_factor), until we have
% less than five prominent peaks in the distribution.

% The goal here is to obtain the modes of the n-modal tempo distributions.
% Next we can use these for various statistics, relative to the tempos of
% the respective songs.
x(x==inf)=[];
x(isnan(x))=[];
if numel(x)<10
    tempo_stats=[];
    return
end

if isempty(varargin)
    plotting_flag = 3;
else
    plotting_flag = varargin{1};
end

[counts,bins]=hist(x,20);
locsm=zeros(10,1);
n_factor=10;
while numel(locsm)>4
    n_factor=n_factor-1;
    %try [~,kpd,dmesh]=kde_botev(x,2^n_factor,0,250);catch;keyboard;end
    [~,kpd,dmesh]=kde_botev(x,2^n_factor,-10,250);
    dmesh=dmesh';
    kpd_rescaled=kpd./sum(kpd);
    switch 2
        case 1
            [pks,locsm]=findpeaks(kpd_rescaled,'MinPeakDistance',3,'MinPeakProminence',.0001);
            pks = pks';
            locsm=dmesh(locsm);
        case 2
            pks = kpd_rescaled(find(diff(diff(kpd_rescaled)>0)==-1)+1);
            locsm = dmesh(find(diff(diff(kpd_rescaled)>0)==-1)+1);
    end
    
    if n_factor==4
        break
    end
    %disp(n_factor)
    %disp(dmesh)
    %disp(bandw)
    
    %subplot(2,1,1)
    %plot(dmesh,kpd_rescaled);hold on;plot(locsm,pks,'s');hold off
    %subplot(2,1,2)
    %plot(dmesh,kpd_rescaled);hold on;plot(dmesh(diff(diff(kpd_rescaled)>0)==-1),kpd_rescaled(diff(diff(kpd_rescaled)>0)==-1),'s');hold off
    %keyboard
end
locsm(pks<.001)=[];
pks(pks<.001)=[];
if isempty(locsm)
    tempo_stats = [];
    return
end
tempo_stats=[locsm pks];
tempo_stats = flipud(sortrows(tempo_stats,2));

kpd_rescaled=kpd_rescaled';
kpd_rescaled=interp(kpd_rescaled,2^9/numel(kpd_rescaled));
dmesh=interp(dmesh,2^9/numel(dmesh));


% Visualize to choose proper parameters for the adaptive kde method.
if plotting_flag==3
    figure(3)
    bar(bins,counts)
    hold on
    plot(dmesh,kpd_rescaled./max(pks)*max(counts),'-r','linewidth',2); %*numel(x)
    plot(locsm,pks./max(pks)*max(counts),'sr'); %*numel(x)
    hold off
    
    % figure(4)
    % plot(dmesh,kpd_rescaled,'-');
    % hold on
    % plot(locsm,pks,'sr');
    % hold off
    % pause
    
    % figure(2)
    % ksdensity(1./cycles(:,2)*60,'bandwidth',7);
    % hold on;
    % kde(1./cycles(:,2)*60,2^8);
    % kde(1./cycles(:,2)*60,2^6);
    % kde(1./cycles(:,2)*60,2^4);
    % hold off
    
    pause
end