function plot_4x8_animated_in_3d(DATA,bodies_labels,wanted_markers,varargin)
%plot_4x8_animated_in_3d
% Video-like plotting.

save_or_plot_flag = 1;

if isempty(varargin)
    sr = 100;
else
    sr = varargin{1};
end

%offset = [14 114 184 26 1 14 48 26];
offset = int16(([19 29 43 16 11 9 16 19]-[7 9 7 1 5 1 1 3])*sr);

ylims = [700 1800;2050 3000;2200 3000;500 1800];
xlims = [0 13500;1000 13000;2000 14500;3000 15000];

win_span_s = 3;

palette = zeros(8,3,4);
for tr=1:8
    for b=1:4
        palette(tr,:,b) = [.3+b/8 1-b/8 .5-b/8]./(1+tr/8);
    end
end
maxt=inf;for tr=1:numel(DATA);maxt = min([maxt size(DATA{tr}.X,1)]);end
maxt = maxt-max(offset);

ds = 1e1/1;
for b = 3:numel(bodies_labels)
    close all
    if save_or_plot_flag == 1
        f = figure('visible','off');
        %set(f,'Position',[1 1 1800 640])
        set(f,'Position',[1 1 1920 1080])
    else
        f = figure;
        %set(f,'Position',[1 1 1800 640])
        set(f,'Position',[1 1 1920 1080])
    end
    
    l=cell(numel(wanted_markers),8);
    p=cell(numel(wanted_markers),8);
    for t = 1:ds:maxt
        for tr = 1:numel(DATA)
            for m = 1:numel(DATA{tr}.col_names)
                for marker = 1:numel(wanted_markers)
                    % Find the indices of the needed markers.
                    if strcmp(DATA{tr}.col_names{m},[bodies_labels{b} wanted_markers{marker}])
                        index = (max([1 t-sr*win_span_s]):t)+offset(tr);
                        x = DATA{tr}.X(index,:,m);
                        if t == 1
                            %plot3(x(:,1)+tr*1.5e3,x(:,2)+2*b*2e3,x(:,3),'-','linewidth',2,'color',palette(tr,:,b))
                            l{marker,tr} = plot3(x(:,1)+tr*1.5e3,x(:,2),x(:,3),'-','linewidth',2,'color',palette(tr,:,b));
                            hold on
                            %plot3(x(end,1)+tr*1.5e3,x(end,2)+2*b*2e3,x(end,3),'o','color','k','MarkerSize',10,'Markerfacecolor',palette(tr,:,b))
                            p{marker,tr} = plot3(x(end,1)+tr*1.5e3,x(end,2),x(end,3),'o','color','k','MarkerSize',10,'Markerfacecolor',palette(tr,:,b));
                            if 1
                                axis off
                                set(gca,'color',[.1 .1 .1])
                                set(gcf,'color','k')
                                %set(gca,'XTickLabels',[])
                                %set(gca,'YTickLabels',[])
                                %set(gca,'ZTickLabels',[])
                            else
                                axis on
                                set(gca,'color','w')
                                set(gcf,'color','w')
                                grid on
                            end
                            
                            xlim(xlims(b,:))
                            ylim(ylims(b,:))
                            zlim([200 1600])
                            
                            set(gca,'view',[30*sawtooth(pi/2+double(t)/sr*2*pi/100,.5) 30+20*sin(double(t)/sr*2*pi/23)])
                            
                        else
                            l{marker,tr}.XData = x(:,1)+tr*1.5e3;
                            l{marker,tr}.YData = x(:,2);
                            l{marker,tr}.ZData = x(:,3);
                            
                            p{marker,tr}.XData = x(end,1)+tr*1.5e3;
                            p{marker,tr}.YData = x(end,2);
                            p{marker,tr}.ZData = x(end,3);

                            %set(gca,'view',[-8 82])
                            %set(gca,'view',[-6 24])
                            %set(gca,'view',[-8.7 27])
                            set(gca,'view',[30*sawtooth(pi/2+double(t)/sr*2*pi/100,.5) 30+20*sin(double(t)/sr*2*pi/23)])
                        end
                    end
                end
            end
            %hold off
        end
        
        if save_or_plot_flag == 1
            set(f, 'InvertHardcopy', 'off')
            print('-dpng','-r100','-loose',['body' num2str(b,'%01.0f') '_frame' num2str(t,'%05.0f') '.png']);
        else
            pause(eps)
        end
        
        fprintf('%8.4f%%',(double(t)/double(maxt))*1e2)
        if mod(t-1,1e2)==0;fprintf('\n');end
    end
end

