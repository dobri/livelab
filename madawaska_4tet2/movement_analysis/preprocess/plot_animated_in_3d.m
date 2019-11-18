function plot_animated_in_3d(DATA)
%plot_animated_in_3d(DATA)
% Video-like plotting.
% Use to quickly get a feel for the data.
% Make sure the ranges are right.
% You can further edit this to make a decent-looking video file.

X = DATA.X;
sf = DATA.sf;

x1=squeeze(X(1,1,:));
y1=squeeze(X(1,2,:));
z1=squeeze(X(1,3,:));
p = plot3(x1,y1,z1,'ko','MarkerSize',5,'MarkerFaceColor','r');
grid on
xlim([-1000 3000])
ylim([-0 3200])
zlim([  500 1500])
set(gca,'view',[-13 74])
for t=2:size(X,1)
    x1=squeeze(X(t,1,:));
    y1=squeeze(X(t,2,:));
    z1=squeeze(X(t,3,:));
    set(p,'XData',x1)
    set(p,'YData',y1)
    set(p,'ZData',z1)
    drawnow
    pause(1/sf)%pause
end
