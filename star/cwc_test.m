function [sig_perc_wtc,mean_rsq,sig_perc_wxy,mean_wxy,anglestrength] = cwc_test(d1,d2,plotting,save_fig2_fname,save_fig3_fname)

wtc_flag = 0;
xwt_flag = 1;

sig_perc_wtc = nan;
mean_rsq = nan;
sig_perc_wxy = nan;
mean_wxy = nan;
anglestrength = nan;

fsz = 30;
dpi = 100;
maxscale = 8;
mc_surr_num = 0;
seriesname={'PP_j' 'PP_i'};

tlim=[min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))];

if plotting==1
    figure(1)
    set(gcf, 'Color', [1 1 1]) %set(gcf,'color',[.5 1 .4])
    subplot(2,1,1);
    wt(d1,'Dj',1/4,'S0',.25,'MaxScale',64,'Mother','morlet','MakeFigure',logical(plotting==1));
    title(seriesname{1});
    set(gca,'xlim',tlim);
    subplot(2,1,2)
    wt(d2,'Dj',1/4,'S0',.25,'MaxScale',64,'Mother','morlet','MakeFigure',logical(plotting==1));
    title(seriesname{2})
    set(gca,'xlim',tlim)
end

if plotting==1
    figure(2)
    set(gcf, 'Color', [1 1 1]) %set(gcf,'color',[.5 .4 1])
end



% XWT
if xwt_flag == 1
    if plotting>0
        if ~isempty(save_fig2_fname)
            figure('Visible', 'off')
            set(gcf, 'Position', [1, 1, dpi*4, dpi*3])
            set(gcf, 'InvertHardcopy', 'off')
        else
            figure(2)
            set(gcf, 'Position', [1+600, 1+400, 300*4, 300*3])
        end
        set(gcf, 'Color', [1 1 1]) %set(gcf, 'Color', [.4 1 .5])
        
        subplot(6,1,5)
        plot(d1(:,1),d1(:,2),'-k','linewidth',2)
        set(gca,'xlim',tlim)
        set(gca,'fontsize',fsz)
        
        subplot(6,1,6)
        plot(d2(:,1),d2(:,2),'-k','linewidth',2)
        xlabel('Time, s')
        set(gca,'xlim',tlim)
        set(gca,'fontsize',fsz)
        
        subplot(6,1,1:4)
    end
    [Wxy,~,scale,coi,sig95] = ...
        xwt(d1,d2,...
        'Dj',1/4,'S0',.25,'MaxScale',maxscale,'Mother','morlet',...
        'MakeFigure',logical(plotting>0));
    if plotting>0
        f=gcf;
        f.Children(1).Position(1)=.8969;
        set(gca,'fontsize',fsz)
        if ~isempty(save_fig3_fname)
            print('-djpeg','-r100',[save_fig2_fname '_' datestr(now,'yyyy-mm-dd-HHMMSS') '.jpeg'])
            close all
        else
            title(['XWT: ' seriesname{1} '-' seriesname{2} ] )
            pause
        end
    end
    
    COI=sig95*0;
    for t=1:size(sig95,2)
        COI(scale<coi(t),t)=1;
    end
    mean_wxy = sum(sum(Wxy.*logical(COI)))/sum(sum(COI));
    sig_perc_wxy = sum(sum((sig95.*logical(COI))>1))./sum(sum(COI));
    [~,anglestrength,~]=anglemean(angle(Wxy.*logical(COI)));
end



% WTC
if wtc_flag == 1
    if plotting>0
        if ~isempty(save_fig3_fname)
            figure('Visible', 'off')
            set(gcf, 'Position', [1, 1, dpi*4, dpi*3])
            set(gcf, 'InvertHardcopy', 'off')
        else
            figure(3)
            set(gcf, 'Position', [1, 1, 300*4, 300*3])
        end
        set(gcf, 'Color', [1 1 1]) %set(gcf, 'Color', [1 .4 .5])
        
        subplot(6,1,5)
        plot(d1(:,1),d1(:,2),'-k','linewidth',2)
        set(gca,'xlim',tlim)
        set(gca,'fontsize',fsz)
        
        subplot(6,1,6)
        plot(d2(:,1),d2(:,2),'-k','linewidth',2)
        set(gca,'xlim',tlim)
        set(gca,'fontsize',fsz)
        xlabel('Time, s')
        
        subplot(6,1,1:4)
    end
    [Rsq,~,scale,coi,sig95] = ...
        wtc(d1,d2,...
        'MonteCarloCount',mc_surr_num,'Dj',1/4,'S0',.25,'MaxScale',maxscale,'Mother','morlet',...
        'MakeFigure',logical(plotting>0));
    
    if plotting>0
        f=gcf;
        f.Children(1).Position(1)=.8969;
        set(gca,'fontsize',fsz)
        if ~isempty(save_fig3_fname)
            print('-djpeg','-r100',[save_fig3_fname '_' datestr(now,'yyyy-mm-dd-HHMMSS') '.jpeg'])
            close all
            %pause
        else
            title(['WTC: ' seriesname{1} '-' seriesname{2} ] )
            pause
        end
    end
    
    COI=sig95*0;
    for t=1:size(sig95,2)
        COI(scale<coi(t),t)=1;
    end
    COI=logical(COI);
    
    sig95=sig95.*COI;
    Rsq=Rsq.*COI;
    
    mean_rsq = sum(sum(Rsq))/sum(sum(COI));
    sig_perc_wtc = sum(sum(sig95>1))./sum(sum(COI));
end
