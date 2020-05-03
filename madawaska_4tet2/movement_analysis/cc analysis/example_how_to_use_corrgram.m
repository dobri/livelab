% This should be self-explanatory. The typical steps in quantifying and
% visualizing cross-correlation between two time series in situations when
% the delay is expected to vary along the length of the recording. Choose
% the range of lags, the window size, and the window overlap based on some
% reasonable task-specific assumptions. Type help corrgram for more detail.
% Look at the author's papers for even more detail.

sr=100;
omega=1.5*pi;
tvec=(0:5e3)'./sr;
x=sin(tvec*omega)+randn(size(tvec))./.5e1;
y=sin(tvec*omega+omega*.2)+randn(size(tvec))./.5e1;
subplot(4,1,1)
plot(tvec,x,tvec,y)

[c,l,t]=corrgram(x,y,round(sr*1),round(sr*5),round(sr*5*.5));
% Or set maxlag based on some externally determined parameter.
% gc_model_order = 1.75; % seconds
% Also, set it to a single window.
% [c,l,t]=corrgram(x,y,round(sr*gc_model_order),numel(x),0);

subplot(4,1,2)
imagesc(t./sr,l./sr,c,[-1 1])
xlabel('Time'), ylabel('Lag, s'), axis xy;
title('Windowed cross-correlation', 'fontweight', 'bold')
colorbar
text(1.2,.5,'C','units','normalized');

[~,ind]=max(c);
cmax=zeros(size(t));
tau_at_cmax=zeros(size(t));
for j=1:numel(t)
    cmax(j)=c(ind(j),j);
    tau_at_cmax(j)=l(ind(j))/sr;
end
subplot(4,1,3)
plot(t./sr,cmax,'-sk')
ylabel('C'), xlabel('Time, s'),
title('The best CC in each window', 'fontweight', 'bold')
subplot(4,1,4)
plot(t./sr,tau_at_cmax,'--sr')
ylabel('\tau, s'), xlabel('Time, s'),
title('The delays that give the best CC in each window', 'fontweight', 'bold')

% Careful with averaging tau!!! 
% The mean of tau doesn't make sense when tau jumps between positive and negative.
fprintf('<C_{max}> = %.3f  at <\\tau> = %.3f s\n',mean(cmax),mean(tau_at_cmax))