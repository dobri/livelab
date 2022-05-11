cd '/home/dobri/logos/livelab/madawaska_4tet2/movement_analysis/data/matlab files/'
load D_Feb142022.mat

dataTrajs={'X_processed'};
sr = 8;
num_trials = 8;
ds = 1;

counter = 0;
for piecei = 1:numel(D)
    for traji = 1:numel(dataTrajs)
        for triali=1:num_trials
            counter = counter + 1;
            figure(1)
            subplot(4,4,counter)
            P = [];
            for pp = 1:4
                x = D{piecei}.(dataTrajs{traji})(pp,:,triali)';
                t = linspace(1,numel(x),numel(x))'./sr;
                [p,f]=pwelch(x,[],[],(.1:.01:4)',sr);
                P = horzcat(P,pow2db(p));
            end
            S = (mean(P.^2,2)-mean(P,2).^2).^.5;
            M = mean(P,2);
            plot(f(2:ds:end),M(2:ds:end),'color',[.8 .1 .1],'linewidth',2)
            set(gca,'color',[.4 .4 .4])
            set(gca,'XTick',0:.5:4)
            text(.9,.9,char(64+counter),'units','normalized')
            grid on
            if counter>12;xlabel('Frequency [Hz]');end
            if any(counter==[1 5 9 13]);ylabel('Power/frequency [dB/Hz]');end
            
            if 0
                color_vec = cool(4);
                figure(2)
                subplot(4,4,counter)
                for pp = 1:4
                    plot3(f(2:ds:end),f(2:ds:end)*0+pp,P(2:ds:end,pp),'color',[color_vec(pp,:) .3],'linewidth',2)
                    hold on
                end
                set(gca,'view',[3 10])
                set(gca,'color',[.5 .5 .5])
                grid on
                hold off
            end
        end
    end
end
