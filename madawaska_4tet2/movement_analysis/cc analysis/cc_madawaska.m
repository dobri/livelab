%% CC_madawaska

% First, load in the data variable 'D'
%load('D.mat')

% Flag for saving data, figs. Set to 1 if you want this loop to save a
% spreadsheet of the data. 0 if no
save_flag=1;
figs_flag=1;
save_wcc_fig=1;
if save_wcc_fig==1
    close all
    figure('Position',[100 100 800*6 800])
    set(gcf,'visible','off')
end

% Specify parameters. 
% Which method are we going to use?
switch 0
    case 0
        method_flag='wcc';
    case 1
        method_flag='cc_and_gcorder';
end

if strcmp(method_flag,'wcc')
    sr8 = 8; % Hz (after downsampling)
    sr100 = 100;
    win_len = 5; % 2.25 seconds
    max_lag = 2; % 1.125 seconds
end
num_trials = numel(D{1}.A);

if strcmp(method_flag,'cc_and_gcorder')
    morders=[D{1}.X_processed_morder,D{1}.X_detrended_processed_morder,...
        D{2}.X_processed_morder,D{2}.X_detrended_processed_morder];
end

%dataTrajs={'X_clean_processed','X_detrended_processed','A'};
%dataTrajs={'X_processed','X_detrended_processed','A'};
dataTrajs={'X_processed','X_detrended_processed'}; %took out A for now



%% CC analysis - position data
counter=0;
for piecei = 1:numel(D)
    for traji = 1:numel(dataTrajs)
        switch dataTrajs{traji}
            case 'A'
                sr = sr100;
            otherwise
                sr = sr8;
        end
        
        % Make vector to store correlation coefficients
        cor_vals = [];
        
        % Make vector to store lags of the max CCs
        max_val_l=[];
        
        % increment the counter
        counter = counter+1;
        
        switch method_flag
            case 'cc_and_gcorder'
                maxlag = morders(counter); %lag is the same as the model order for gc
                window=length(D{piecei}.(dataTrajs{traji})); %Whole length of the piece.
                overlap=0; %No overlap
            case 'wcc'
                maxlag=round(max_lag*sr); % 2 seconds
                window=round(win_len*sr); % 5 seconds
                overlap=round(window/2); % half a window overlap
        end
        
        for triali=1:num_trials
            switch 1 
                case 0
                    %salient_points = (30:30:max(t/sr))';
                    salient_points = cumsum(randi(10,10,1)+30-5);
                case 1
                    % Read a clicktrack from a text file such as the marked labels exported from Audacity.
                    S = readtable(['p' num2str(piecei) 't' num2str(triali) '.txt']);
                    salient_points = S.Var1;
                    salient_points_markers=S.Var3;
            end
            
            %take the maximum unsigned CC coefficient for each of the 6 possible
            %pairs of musicians for each trial
            switch method_flag
                case 'cc_and_gcorder'
                    % lag is the same as the model order for gc
                    % single trial-long window
                    mcounter = 0;
                    for row=1:4
                        for col=1:4
                            if col>row
                                mcounter = mcounter + 1;
                                [c,l]=xcov(D{piecei}.(dataTrajs{traji})(row,:,triali),D{piecei}.(dataTrajs{traji})(col,:,triali),morders(counter),'coef'); %'coef' - normalizes the sequence
                                cor_vals(mcounter+6*(triali-1),1,piecei)=max(abs(c));
                            end
                        end
                    end
                case 'wcc'
                    % max(A,[],'all') works from R2018b after.
                    fcounter = 0;
                    for row=1:4
                        for col=1:4
                            if col>row
                                fcounter = fcounter + 1;
                                if traji==3
                                    x = D{piecei}.A{triali}(row,:)';
                                    y = D{piecei}.A{triali}(col,:)';
                                else
                                    x = D{piecei}.(dataTrajs{traji})(row,:,triali);
                                    y = D{piecei}.(dataTrajs{traji})(col,:,triali);
                                end
                                [wcc,l,t]=corrgram(x,y,maxlag,window,overlap);
                                cor_vals(fcounter+6*(triali-1),1,piecei)=max(max(abs(wcc)));
                                [indexR,indexC]=find(abs(wcc)==max(max(abs(wcc)))); %find the lag corresponding to the wcc value
                                max_val_l(fcounter+6*(triali-1),1,piecei)=indexR;
                                
                                if figs_flag
                                    % This makes a pseudo-time vector that
                                    % runs in units of time defined by the
                                    % time points in the salient_points vector.
                                    tn = time_interpolate(t./sr,salient_points);
                                    stn = time_interpolate(salient_points,salient_points);
                                    
                                    %corrgram(D{piecei}.(dataTrajs{traji})(row,:,triali),D{piecei}.(dataTrajs{traji})(col,:,triali),maxlag,window,overlap)
                                    if fcounter == 1 %I need to do this or else no plot appears
                                        figure
                                    end 
                                    
                                    subplot(num_trials,6,fcounter+6*(triali-1)) %note to self - change so there is a separate plot for each pair
                                    %imagesc(tn,l,flip(wcc),[-1 1]);colorbar
                                    imagesc(t./sr,l,flip(wcc),[-1 1]);colorbar
                                    xtickangle(30)
                                    colorbar off
                                    %if fcounter~=6
                                    %    colorbar off
                                    %end
                                    if fcounter==1
                                        lags_show = unique(sort([0;l(1);l(end);l(1:(sr/4):end)]));
                                        set(gca,'YTick',lags_show)
                                        set(gca,'YTickLabel',flipud(lags_show./sr))
                                        ylabel('Lag, s')
                                    else
                                        set(gca,'YTick',[])
                                        ylabel('')
                                    end
                                    %set(gca,'XTick',1:numel(salient_points))
                                    set(gca,'XTick',salient_points)
                                    set(gca,'XTickLabel',round(salient_points,2))
                                    line([salient_points';salient_points'],[salient_points';salient_points']*0+[min(lags_show);max(lags_show)],'color','r','linestyle','--')
                                    %line([stn';stn'],[stn';stn']*0+[min(lags_show);max(lags_show)],'color','r')
                                    %set(gca,'XTick',t(1:10:end))
                                    %set(gca,'XTickLabel',round(t(1:10:end)./sr))
                                    %if triali==size(D{piecei}.(dataTrajs{traji}),3)
                                    %    xlabel('Time [seconds or sections]')
                                    %end
                                    if triali==1
                                        if row==1
                                            row_lab='cello';
                                        elseif row==2
                                            row_lab='viola';
                                        elseif row==3
                                            row_lab='violin1';
                                        else
                                            row_lab='violin2';
                                        end
                                        
                                        if col==1
                                            col_lab='cello';
                                        elseif col==2
                                            col_lab='viola';
                                        elseif col==3
                                            col_lab='violin1';
                                        else
                                            col_lab='violin2';
                                        end
                                        %title(['Pair ' num2str(row) '-' num2str(col)])
                                        title([row_lab '-' col_lab])

                                    else
                                        title('')
                                    end
                                end
                            end
                        end
                    end
            end
        end
        
        if figs_flag
            if save_wcc_fig == 1 % print to file
                set(gca(), 'LooseInset', get(gca(), 'TightInset'));
                %print(gcf,'-dpng','-r300','-loose',['wcc_' 'score' num2str(piecei) '_dim' num2str(traji) ...
                %    '_tr' num2str(triali) '_' datestr(now,'yymmdd-HHMMSS') '.png']);
                print(gcf,'-dpng','-r100','-loose',['wcc_' 'score' num2str(piecei) '_dim' num2str(traji) ...
                    '_' datestr(now,'yymmdd-HHMMSS') '.png']);
            else
                pause % just inspect on the screen
            end
        end
        
        label_cc=[dataTrajs{traji},'_',method_flag];
        label_cc_lags=[dataTrajs{traji},'_',method_flag,'_lag'];

        D{piecei}.(label_cc)=cor_vals(:,:,piecei);
        D{piecei}.(label_cc_lags)=max_val_l(:,:,piecei);

    end
end

if save_wcc_fig == 1
    set(gcf,'visible','on')
    close all
end


figure
for p=1:2
    subplot(2,2,1+(p-1)*2)
    boxplot(D{p}.X_processed_wcc,reshape(meshgrid(1:8,1:6),1,[])')
    subplot(2,2,2+(p-1)*2)
    boxplot(D{p}.X_detrended_processed_wcc,reshape(meshgrid(1:8,1:6),1,[])')
end


figure
for p=1:2
    subplot(2,2,1+(p-1)*2)
    boxplot(D{p}.A_wcc,reshape(meshgrid(1:8,1:6),1,[])')
end


%% Save data
if save_flag==1
    
    %TURN THIS INTO A trajectory LOOP once I FIX the piece 2 labels
    
    %Position
    
    %Reconfigure cor_vals
    cor_vals_reconfig=[D{1}.X_processed_wcc;D{2}.X_processed_wcc]; %fix this once I get piece 2 markers
    
    %Reconfigure lags
    cor_lags_reconfig=[D{1}.X_processed_wcc_lag;D{2}.X_processed_wcc_lag]; %fix this once I get piece 2 markers
    
    %Make vector for pair
    pair=repmat([1:6]',16,1);
    
    %Make vector with the pair names
    pair_names={'cello-viola','cello-v1','cello-v2','viola-v1','viola-v2','v1-v2'}';
    pair_names=[pair_names;pair_names;pair_names;pair_names;pair_names;pair_names;pair_names;pair_names;pair_names;pair_names;pair_names;pair_names;pair_names;pair_names;pair_names;pair_names];
    
    %Make vector for piece
    piece=repelem([1;2],48);
    
    %Make vector for trial
    trial=[1+zeros(6,1);2+zeros(6,1);3+zeros(6,1);4+zeros(6,1);...
        5+zeros(6,1);6+zeros(6,1);7+zeros(6,1);8+zeros(6,1);];
    trial=repmat(trial,2,1);
    
    %Make vector for trial_collapsed
    trial_collapsed = repmat(repelem([1:4]',12),2,1);
    
    
    %Make vector for condition
    condition1=[1+zeros(6,1);2+zeros(6,1);1+zeros(6,1);2+zeros(6,1);...
        1+zeros(6,1);2+zeros(6,1);1+zeros(6,1);2+zeros(6,1)];
    condition2=flip(condition1);
    condition=[condition1;condition2];
    
    
    T=table(pair, pair_names, cor_vals_reconfig,cor_lags_reconfig, condition,trial,trial_collapsed,piece);
    filename='mada_cc_position_window.xlsx';
    T.Properties.VariableNames = {'pair','pair_names', 'wcc', 'lag','condition','trial', 'trial_collapsed','piece'};
    writetable(T,filename);
        
    %Detrended
    %Reconfigure cor_vals
    cor_vals_reconfig_detrended=[D{1}.X_detrended_processed_wcc;D{2}.X_detrended_processed_wcc]; 
    T2=table(pair, pair_names, cor_vals_reconfig_detrended,condition,trial,trial_collapsed,piece);
    filename2='mada_cc_position_detrended_window.xlsx';
    T2.Properties.VariableNames = {'pair','pair_names', 'wcc', 'condition','trial', 'trial_collapsed','piece'};
    writetable(T2,filename2);
    
    %Acceleration
    cor_vals_reconfig_A=[D{1}.A_wcc;D{2}.A_wcc];
    T3=table(pair, pair_names, cor_vals_reconfig_A,condition,trial,trial_collapsed,piece);
    filename3='mada_cc_A_window.xlsx';
    T3.Properties.VariableNames = {'pair','pair_names', 'wcc', 'condition','trial', 'trial_collapsed','piece'};
    writetable(T3,filename3);
    

end

