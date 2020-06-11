function sig_perc = cwc_test(d1,d2,plotting,save_fig3_fname)
%cd ~
%addpath('~/logos/mcmc/livelab/star/wavelet-coherence/')

mc_surr_num = 2e2;
seriesname={'PP_j' 'PP_i'};

%ds=4;
%fr = lvmeta.fr/ds; %100
%d1=[lvtime{1}(1:ds:end) lvposu{1}(1:ds:end,1,1)];
%d2=[lvtime{1}(1:ds:end) lvposu{1}(1:ds:end,2,1)];

if plotting==1
    figure(1)
    set(gcf,'color',[.5 1 .4])
end
tlim=[min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))];

if plotting==1
    subplot(2,1,1);
end
wt(d1,'Dj',1/4,'S0',.25,'MaxScale',64,'Mother','morlet','MakeFigure',logical(plotting==1));
if plotting==1
    title(seriesname{1});
    set(gca,'xlim',tlim);
end

if plotting==1
    subplot(2,1,2)
end
wt(d2,'Dj',1/4,'S0',.25,'MaxScale',64,'Mother','morlet','MakeFigure',logical(plotting==1));
if plotting==1
    title(seriesname{2})
    set(gca,'xlim',tlim)
end

if plotting==1
    figure(2)
    set(gcf,'color',[.5 .4 1])
end
xwt(d1,d2,'Dj',1/4,'S0',.25,'MaxScale',32,'Mother','morlet','MakeFigure',logical(plotting==1));
if plotting==1
    title(['XWT: ' seriesname{1} '-' seriesname{2} ] )
end

% wtc(d1,d2,'MonteCarloCount',1e4,'Dj',1/4,'S0',2,'MaxScale',20,'Mother','morlet')
%{
  .         Dj: Octaves per scale (default: '1/12')
  .         S0: Minimum scale
  .         J1: Total number of scales
  .         Mother: Mother wavelet (default 'morlet')
  .         MaxScale: An easier way of specifying J1
%}
if plotting>0
    if ~isempty(save_fig3_fname)
        figure('Visible', 'off')
    else
        figure(3)
    end
    set(gcf, 'InvertHardcopy', 'off')
    set(gcf, 'Color', [1 .4 .5])
    set(gcf, 'Position', [1, 1, 300*4, 300*3])
    
    subplot(6,1,5)
    plot(d1(:,1),d1(:,2),'-k','linewidth',2)
    subplot(6,1,6)
    plot(d2(:,1),d2(:,2),'-k','linewidth',2)
    xlabel('Time, s')
    
    subplot(6,1,1:4)
end
[Rsq,period,scale,coi,sig95] = ...
    wtc(d1,d2,...
    'MonteCarloCount',mc_surr_num,'Dj',1/4,'S0',.25,'MaxScale',32,'Mother','morlet',...
    'MakeFigure',logical(plotting>0));
if plotting>0
    if ~isempty(save_fig3_fname)
        print('-djpeg','-r100',[save_fig3_fname '_' datestr(now,'yyyy-mm-dd-HHMMSS') '.jpeg'])
        close all
    else
        title(['WTC: ' seriesname{1} '-' seriesname{2} ] )
    end
end

% imagesc(d1(:,1),period,sig95-1>1)
sig_perc = sum(sum(sig95>1))./numel(sig95);
