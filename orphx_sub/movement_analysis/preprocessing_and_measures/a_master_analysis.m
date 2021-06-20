clear

% Move to where the data is stored, start reading-in parameters.
% i.e. '/media/dobri/disk2/mcmc_large_data/orphx'
data_path = input('Me: Where is my money?\nYou: ');
cd(data_path)
load(fullfile(data_path,'orphx_set2_preprocessed_and_cleaned_2021-06-19.mat'))

% Set this to visually inspect the multiple stages of data processing.
plotting_flag = 0;

% What is the frequency range of head-bobbing? First I used .2 to 4 Hz.
[b,a]=butter(2,[.1/(.802*ll.sf/2) 4/(.802*ll.sf/2)]);

mid_processing.pca.eigen = cell(ll.num_trial,ll.num_pp);
mid_processing.pca.r2 = nan(ll.num_trial,3,ll.num_pp);
mid_processing.pca.r2_pc2_over_pc1 = nan(ll.num_trial,1,ll.num_pp);
mid_processing.pca.pc = cell(ll.num_trial,ll.num_pp);
mid_processing.cycl = cell(ll.num_trial,ll.num_pp);

% Return with Trial, PPNUM, PPID, CONDITION, A, BPM, H, R2PC2/R2PC1, ...
rez.meta = nan(ll.num_trial*ll.num_pp,4);
rez.amp = nan(ll.num_trial*ll.num_pp,1);
rez.tempo_mode = nan(ll.num_trial*ll.num_pp,1);
rez.pc2_over_pc1 = nan(ll.num_trial*ll.num_pp,1);
rez.h = nan(ll.num_trial*ll.num_pp,1);
rez.q1 = nan(ll.num_trial*ll.num_pp,1);
rez.q2 = nan(ll.num_trial*ll.num_pp,1);
rez.v = nan(ll.num_trial*ll.num_pp,1);
rez.pl = nan(ll.num_trial*ll.num_pp,1);
rez.osc_prop = nan(ll.num_trial*ll.num_pp,1);

% Analyze per trial.
% Most of this analysis ends up being irrelevant. It's heavily focused on
% oscillations. It seems like in the present experiment participant didn't
% respond strongly and directly to the beat, hence their body sway was not
% oscillating as it would if they were dancing. That's why a simple path
% length analysis looks to make better sense. At a certain point, it could
% be useful to explicitly make that point. In the STAR study, participants'
% postural control policy was that of a pendulum. Here, brownian particle.
% Ask me for the relevant literature if necessary.
for tr=1:ll.num_trial
    for pp=1:ll.num_pp
        table_row = (tr-1)*ll.num_pp+pp;
        fprintf('%4.0f,%4.0f\n',tr,pp)
        rez.meta(table_row,:) = [tr pp ll.IDSnum(pp) abs(mod(tr,2)-1)+1];
        D=ll.D{tr}(:,:,pp);

        % Skip if the nan prop > .5.
        if sum(isnan(D(:,1)))/size(D,1)>.5
            continue
        end
        
        % Crop the nans. It seems like there isn't a great number of
        % intermittent nan-intervals within a participant's trial. Hence, 
        % don't worry about difficult stitching procedures. Just crop nans.
        D(sum(isnan(D),2)>0,:)=[];
        time=((1:size(D,1))')./ll.sf;
        
        % Filter & PCA.
        [eigenvec,temp,~,~,explained]=pca(filtfilt(b,a,D));
        mid_processing.pca.eigen{tr,pp} = eigenvec;
        mid_processing.pca.pcx{tr,pp} = temp;
        mid_processing.pca.r2(tr,:,pp) = explained;
        mid_processing.pca.r2_pc2_over_pc1(tr,:,pp) = explained(2)/explained(1);
        rez.pc2_over_pc1(table_row) = mid_processing.pca.r2_pc2_over_pc1(tr,1,pp);

        if plotting_flag == 1
            figure(1)
            subplot(2,1,1)
            plot(time,D)
            hold on
            plot(time,filtfilt(b,a,D))
            hold off
            title('Movement, raw and filtered')
            
            subplot(2,1,2)
            plot(time,temp)
            title('PCs')
            
            pause
        end
        
        % Oscillations analysis: Peak times, periods, and amplitudes.
        % The oscillation energy H (sum across both PCs).
        mid_processing.cycl{tr,pp} = cycles_amp_freq_energy(temp(:,1:2),time,ll.sf,plotting_flag*2);
        switch 0 % for debugging, observe the raw figures.
            case 2
                print('-djpeg','-r72',['raw_mov_amp_h_etc_TR' num2str(tr,'%02.0f') '_PP' num2str(pp,'%2.0f') '_' datestr(now,'yyyy-mm-dd-HHMMSS') '.jpeg'])
        end

        % Summarize the mov variables as means or medians.
        if ~isempty(mid_processing.cycl{tr,pp})
            rez.osc_prop(table_row,1) = mid_processing.cycl{tr,pp}.osc_prop;
            rez.q1(table_row) = mean(mid_processing.cycl{tr,pp}.q(300:end,1));
            rez.q2(table_row) = mean(mid_processing.cycl{tr,pp}.q(300:end,2));
            rez.amp(table_row) = median(mid_processing.cycl{tr,pp}.amp_interp(300:end,1));
            rez.h(table_row) = log(mean(mid_processing.cycl{tr,pp}.h(300:end,end)));
            %rez.h(table_row) = median(mid_processing.cycl{tr,pp}.h(300:end,end));
            
            mid_processing.cycl{tr,pp}.v = dot(diff(D(300:end,:)),diff(D(300:end,:)),2).^.5./(1/ll.sf);
            rez.v(table_row,1) = median(mid_processing.cycl{tr,pp}.v);
            rez.pl(table_row,1) = sum(dot(diff(D(300:end,:)),diff(D(300:end,:)),2).^.5)/(size(D,1)/ll.sf/60);
        end

        % BPM mode, as the location of the largest peak in the kde.
        % Pool PC1 and 2 together. If needed, the analysis can be split.
        switch 1
            case 1
                bpms = [60*mid_processing.cycl{tr,pp}.cycles{1}(:,2); 60*mid_processing.cycl{tr,pp}.cycles{2}(:,2)];
                mid_processing.cycl{tr,pp}.tempo_stats = tempo_modes_with_adaptive_kde(bpms,plotting_flag*3);
                if ~isempty(mid_processing.cycl{tr,pp}.tempo_stats)
                    rez.tempo_mode(table_row) = mid_processing.cycl{tr,pp}.tempo_stats(1);
                end
            case 2
                for d=1:2
                    mid_processing.cycl{tr,pp}.tempo_stats{d} = tempo_modes_with_adaptive_kde(60*mid_processing.cycl{tr,pp}.cycles{d}(:,2),plotting_flag*3);
                end
        end
    end
end

rez.out = array2table([rez.meta rez.tempo_mode rez.pc2_over_pc1 rez.amp rez.h rez.q1 rez.q2 rez.v rez.pl rez.osc_prop]);
rez.out.Properties.VariableNames = {'Trial','pp','ID','Condition','Amp','R2PC2byPC1','BPMmode','H','q1','q2','Speed','PathLength','OscTimeProp'};
if plotting_flag
    close all
    for dv=5:13
        subplot(3,3,dv-4)
        boxplot(table2array(rez.out(:,dv)),table2array(rez.out(:,4)))
        title(rez.all.labels{dv})
    end
end

%{
save(fullfile(data_path,['rez_' datestr(now,'yyyy-mm-dd') '.mat']),'rez')
writetable(rez.out,fullfile(data_path,['rez_' datestr(now,'yyyy-mm-dd') '.csv']))
%}
