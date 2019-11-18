plotting = 0;

% V Manually mark the points in time when participants change spatial location.
% V Zero-center the separate periods.
% V Broken into trials. 
%cd('/home/dobri/mcmc/orphx/manalysis')
%addpath('~/mcmc/orphx/manalysis/')
load('orphx_set2_preprocessed_and_cleaned_2019-08-29.mat')

% What is the frequency range of head-bobbing? First I used .2 to 4 Hz.
[b,a]=butter(2,[.1/(.802*ll.sf/2) 4/(.802*ll.sf/2)]);

mid_processing.pca.eigen = cell(18,43);
mid_processing.pca.r2 = nan(18,3,43);
mid_processing.pca.r2_pc2_over_pc1 = nan(18,1,43);
mid_processing.pca.pc = cell(18,43);
mid_processing.cycl = cell(18,43);

% Return with Trial, PPNUM, PPID, CONDITION, A, BPM, H, R2PC2/R2PC1, ...
rez.meta = nan(18*43,4);
rez.amp = nan(18*43,1);
rez.tempo_mode = nan(18*43,1);
rez.pc2_over_pc1 = nan(18*43,1);
rez.h = nan(18*43,1);
rez.q1 = nan(18*43,1);
rez.q2 = nan(18*43,1);
rez.v = nan(18*43,1);
rez.pl = nan(18*43,1);
rez.osc_prop = nan(18*43,1);

% V Analyze per trial. tr=1;pp=1;
for tr=1:18
    %discarded=[];
    for pp=1:43
        fprintf('%4.0f,%4.0f\n',tr,pp)
        rez.meta(((tr-1)*43+pp),:) = [tr pp ll.IDSnum(pp) abs(mod(tr,2)-1)+1];
        D=ll.D{tr}(:,:,pp);
        % disp(size(D))
        % V Skip if the nan prop > .5.
        if sum(isnan(D(:,1)))/size(D,1)>.5
            continue
        end
        % V Crop the nans. (Preserve the nan index?)
        D(sum(isnan(D),2)>0,:)=[];
        time=((1:size(D,1))')./ll.sf;
        

        % It seems like there isn't a great number of intermittent
        % nan-intervals within a participant's trial. Hence, don't worry
        % about a difficult stitching and cleaning procedure. Just crop nans.
        % Of course, remove participants with a lot of data missing (>50%).
        %         if sum(isnan(D(:,1)))/size(D,1)>.5
        %             discarded=[discarded,pp];
        %             continue
        %         end
        %         if sum(isnan(D(:,1)))/size(D,1)>.01
        %             subplot(2,1,1)
        %             plot(D)
        %             hold on
        %             plot((sum(isnan(D),2)>0)*100)
        %             hold off
        %             D(sum(isnan(D),2)>0,:)=[];
        %             subplot(2,1,2)
        %             plot(D)
        %             pause
        %         end

        
        % Filter & PCA.
        [eigenvec,temp,~,~,explained]=pca(filtfilt(b,a,D));
        mid_processing.pca.eigen{tr,pp} = eigenvec;
        mid_processing.pca.pcx{tr,pp} = temp;
        mid_processing.pca.r2(tr,:,pp) = explained;
        mid_processing.pca.r2_pc2_over_pc1(tr,:,pp) = explained(2)/explained(1);
        rez.pc2_over_pc1((tr-1)*43+pp) = mid_processing.pca.r2_pc2_over_pc1(tr,1,pp);

        if plotting==1
            figure(1)
            subplot(2,1,1)
            plot(time,D)
            hold on
            plot(time,filtfilt(b,a,D))
            hold off
            subplot(2,1,2)
            plot(time,temp)
            pause
        end
        
        % V Oscillations analysis: Peak times, periods, and amplitudes.
        % V The oscillation energy H (sum across both PCs).
        % cycl.cycles, cycl.peak_loc, cycl.peak_amp, cycl.iois,
        % cycl.amp_interp, cycl.freq_interp, cycl.h, cycl.q
        mid_processing.cycl{tr,pp} = cycles_amp_freq_energy(temp(:,1:2),time,ll.sf,2);
        switch 2
            case 2
                print('-djpeg','-r72',['raw_mov_amp_h_etc_TR' num2str(tr,'%02.0f') '_PP' num2str(pp,'%2.0f') '_' datestr(now,'yyyy-mm-dd-HHMMSS') '.jpeg'])
        end

        % ! Insert an image-printing option for the raw cycle analysis.

        if ~isempty(mid_processing.cycl{tr,pp})
            rez.q1((tr-1)*43+pp) = mean(mid_processing.cycl{tr,pp}.q(300:end,1));
            rez.q2((tr-1)*43+pp) = mean(mid_processing.cycl{tr,pp}.q(300:end,2));
            rez.h((tr-1)*43+pp) = log(mean(mid_processing.cycl{tr,pp}.h(300:end,end)));
            rez.amp((tr-1)*43+pp) = median(mid_processing.cycl{tr,pp}.amp_interp(300:end,1));
            
            mid_processing.cycl{tr,pp}.v = dot(diff(D(300:end,:)),diff(D(300:end,:)),2).^.5./(1/ll.sf);
            rez.v((tr-1)*43+pp,1) = median(mid_processing.cycl{tr,pp}.v);
            rez.osc_prop((tr-1)*43+pp,1) = mid_processing.cycl{tr,pp}.osc_prop;
            rez.pl((tr-1)*43+pp,1) = sum(dot(diff(D(300:end,:)),diff(D(300:end,:)),2).^.5)/(size(D,1)/ll.sf/60);
        end

        
        % V BPM mode. Pool pc1 and 2 together. If needed, the analysis can be split.
        % V BPM as the location of the largest peak in the kde, from among both PC1 and PC2.
        switch 1
            case 1
                bpms = [60*mid_processing.cycl{tr,pp}.cycles{1}(:,2); 60*mid_processing.cycl{tr,pp}.cycles{2}(:,2)];
                mid_processing.cycl{tr,pp}.tempo_stats = tempo_modes_with_adaptive_kde(bpms,0);
                if ~isempty(mid_processing.cycl{tr,pp}.tempo_stats)
                    rez.tempo_mode((tr-1)*43+pp) = mid_processing.cycl{tr,pp}.tempo_stats(1);
                end
            case 2
                for d=1:2
                    subplot(1,2,d)
                    mid_processing.cycl{tr,pp}.tempo_stats{d} = tempo_modes_with_adaptive_kde(60*mid_processing.cycl{tr,pp}.cycles{d}(:,2),0);
                end
        end
        % Later, if it's necessary to show the bpm modes across trials and
        % participants, find the following among the star scripts.
        % fig_kde_prob_densities_tempo_nice(kde_tempos)
        %print('-djpeg','-r600',['tempo_distr_kde_' datestr(now,'yyyy-mm-dd-HHMMSS') '.jpeg'])
        %print('-depsc',['tempo_distr_kde_' datestr(now,'yyyy-mm-dd-HHMMSS') '.eps'])
    end
    %    discarded
end

rez.all.labels = {'Trial','pp','ID','Cond','Amp','R2PC2byPC1','BPM-mode','H','q1','q2','Speed','PathLength','OscTimeProp'};
rez.all.data = [rez.meta rez.tempo_mode rez.pc2_over_pc1 rez.amp rez.h rez.q1 rez.q2 rez.v rez.pl rez.osc_prop];
for dv=5:13;subplot(3,3,dv-4);boxplot(rez.all.data(:,dv),rez.all.data(:,4));title(rez.all.labels{dv});end

%{
print('-djpeg','-r300',['rez_dvs_' datestr(now,'yyyy-mm-dd-HHMMSS') '.jpeg'])

save(fullfile(pwd,['rez_' datestr(now,'yyyy-mm-dd-HHMMSS') '.mat']),'rez')

fid=fopen(fullfile(pwd,['rez_' datestr(now,'yyyy-mm-dd-HHMMSS') '.csv']),'a');
for l=1:numel(rez.all.labels);
fprintf(fid,'%s,',rez.all.labels{l});
end;
fprintf(fid,'\n');
for k=1:size(rez.all.data,1);
fprintf(fid,'%12.4f,',rez.all.data(k,:));
fprintf(fid,'\n');
end;
fclose(fid);
%}
