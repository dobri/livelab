function plot_4x8_animated_in_3d(DATA,bodies_labels,wanted_head_markers,varargin)
%plot_4x8_animated_in_3d
% Video-like plotting.

offset = [14 114 184 26 1 14 48 26];

if isempty(varargin)
    sf = 100;
else
    sf = varargin{1};
end

for tr=1:8
    for b=1:4
        palette(tr,:,b) = [.3+b/8 1-b/8 .5-b/8]./(1+tr/8);
    end
end
maxt=inf;for tr=1:numel(DATA);maxt = min([maxt size(DATA{tr}.X,1)]);end
maxt = maxt-max(offset);

for t = 1:1e1:maxt
    for tr = 1:numel(DATA)
        for b = 1:numel(bodies_labels)
            % Find the indices of the needed markers.
            for m = 1:numel(DATA{tr}.col_names)
                for marker = 1:numel(wanted_head_markers)
                    if strcmp(DATA{tr}.col_names{m},[bodies_labels{b} wanted_head_markers{marker}])
                        index = (max([1 t-sf*4]):t)+offset(tr);
                        x = DATA{tr}.X(index,:,m);
                        plot3(x(:,1)+tr*1.5e3,x(:,2)+2*b*2e3,x(:,3),'-','linewidth',2,'color',palette(tr,:,b))
                        hold on
                        plot3(x(end,1)+tr*1.5e3,x(end,2)+2*b*2e3,x(end,3),'o','color','k','MarkerSize',10,'Markerfacecolor',palette(tr,:,b))
                    end
                end
            end
        end
    end
    hold off
    set(gcf,'color','k')
    set(gca,'color',[.1 .1 .1])
    grid on
    xlim([-1000 15000])
    ylim([  500 20000])
    zlim([  200 1500])
    set(gca,'view',[-8 82])
    %pause(eps)
    set(gcf, 'InvertHardcopy', 'off')
    print(gcf,'-dpng','-r100','-loose',['frame' num2str(t,'%05.0f') '.png']);
end

