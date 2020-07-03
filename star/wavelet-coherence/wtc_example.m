cd ~/logos/mcmc/livelab/star/wavelet-coherence/

seriesname={'AO' 'BMI'};
d1=load('faq/jao.txt');
d2=load('faq/jbaltic.txt');

d2(:,2)=boxpdf(d2(:,2));

figure('color',[1 1 1])
tlim=[min(d1(1,1),d2(1,1)) max(d1(end,1),d2(end,1))];
subplot(2,1,1);
wt(d1);
title(seriesname{1});
set(gca,'xlim',tlim);
subplot(2,1,2)
wt(d2)
title(seriesname{2})
set(gca,'xlim',tlim)

figure('color',[1 1 1])
xwt(d1,d2)
title(['XWT: ' seriesname{1} '-' seriesname{2} ] )

figure('color',[1 .4 .5])
% wtc(d1,d2,'MonteCarloCount',1e4,'Dj',1/4,'S0',2,'MaxScale',20,'Mother','morlet')
%{
  .         Dj: Octaves per scale (default: '1/12')
  .         S0: Minimum scale
  .         J1: Total number of scales
  .         Mother: Mother wavelet (default 'morlet')
  .         MaxScale: An easier way of specifying J1
%}
[Rsq,period,scale,coi,sig95] = ...
    wtc(d1,[d2(:,1) (d2(:,2))],'MonteCarloCount',1e3,'Dj',1/4,'S0',2,'MaxScale',64,'Mother','morlet','MakeFigure',true);
title(['WTC: ' seriesname{1} '-' seriesname{2} ] )

imagesc(d1(:,1),(1:(1/4):6),((sig95-1)>.05))
