function plot_animated_in_3d(X,varargin)
%plot_animated_in_3d(X,sf)
% Video-like plotting.
% X is time by dim by marker.
% Use to quickly get a feel for the data.
% Make sure the ranges are right.
% You can further edit this to make a decent-looking video file.

if isempty(varargin)
    sf = 100;
else
    sf = varargin{1};
end

if isempty(varargin)
    real_units_time_plot_flag = 0;
else
    real_units_time_plot_flag = varargin{2};
end

if numel(varargin)<3
    skip_seconds = 0;
else
    skip_seconds = varargin{3};
end

if real_units_time_plot_flag == 1
    ds = 1;
else
    ds=10;
end

x1=squeeze(X(1,1,:));
y1=squeeze(X(1,2,:));
z1=squeeze(X(1,3,:));
p = plot3(x1,y1,z1,'ko','MarkerSize',5,'MarkerFaceColor','r');
tx = text(0,1,1,num2str(0,'%8.2f secs'),'units','normalized','fontsize',10);
grid on
xlim([-1000 3000])
ylim([-0 3200])
zlim([  500 1500])
set(gca,'view',[-13 74])
for t=(skip_seconds*sf+2):ds:size(X,1)
    tic
    x1=squeeze(X(t,1,:));
    y1=squeeze(X(t,2,:));
    z1=squeeze(X(t,3,:));
    set(p,'XData',x1)
    set(p,'YData',y1)
    set(p,'ZData',z1)
    set(tx,'String',num2str(t/sf,'%8.2f secs'))
    drawnow
    if real_units_time_plot_flag
        pause(1/sf-toc)
    else
        pause(eps)
    end
end
