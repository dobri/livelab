%% CC_madawaska

% First, load in the data variable 'D'

% Load data
%load('D.mat')

% Flag for saving data. Set to 1 if you want this loop to save a
% spreadsheet of the data. 0 if no

save_flag=0;
figs_flag=1;
save_wcc_fig=1;
if save_wcc_fig==1
    close all
    figure('Position',[100 100 1600 800])
    set(gcf,'visible','off')
end

%Specify parameters
%Which method are we going to use?
switch 0
    case 0
        method_flag='wcc';
    case 1
        method_flag='cc_and_gcorder';
end

if strcmp(method_flag,'wcc')
    sr = 8; % Hz. After downsampling?
    win_len = 5; % seconds
    max_lag = 2; % seconds
end

if strcmp(method_flag,'cc_and_gcorder')
    morders=[D{1}.X_processed_morder,D{1}.X_detrended_processed_morder,...
    D{2}.X_processed_morder,D{2}.X_detrended_processed_morder];
end

dataTrajs={'X_processed','X_detrended_processed'};

%% CC analysis - position data
counter=0;
for piecei = 1:numel(D)
    for traji = 1:numel(dataTrajs)
        % Preallocate vector to store correlation coefficients
        cor_vals=zeros(6*size(D{piecei}.(dataTrajs{traji}),3),1,numel(D));
        
        % Specify the same parameters that Andrew used
        counter=counter+1;
        
        switch method_flag
            case 'cc_and_gcorder'
                maxlag = morders(counter); %lag is the same as the model order for gc
                window=length(D{piecei}.(dataTrajs{traji})); %Whole length of the piece. 
                overlap=0; %No overlap
            case 'wcc'
                maxlag=round(max_lag*sr); % 2 seconds
                window=round(win_len*sr); % 5 seconds 
                overlap=round(window/5); % half a window overlap
        end

        for triali=1:size(D{piecei}.(dataTrajs{traji}),3)
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
                                [c,l]=xcov(D{piecei}.(dataTrajs{traji})(row,:,triali),D{piecei}.(dataTrajs{traji})(col,:,triali),morders(counter),'coef'); %only thing I'm not sure about is 'coef' - normalizes the sequence 
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
                                [wcc,l,t]=corrgram(D{piecei}.(dataTrajs{traji})(row,:,triali),D{piecei}.(dataTrajs{traji})(col,:,triali),maxlag,window,overlap);
                                cor_vals(fcounter+6*(triali-1),1,piecei)=max(max(abs(wcc)));
                                if figs_flag
                                    subplot(2,3,fcounter)
                                    corrgram(D{piecei}.(dataTrajs{traji})(row,:,triali),D{piecei}.(dataTrajs{traji})(col,:,triali),maxlag,window,overlap)
                                    xtickangle(30)
                                    set(gca,'YTick',l(1:4:end))
                                    set(gca,'YTickLabel',l(1:4:end)./sr)
                                    ylabel('Lag, s')
                                    set(gca,'XTick',t(1:10:end))
                                    set(gca,'XTickLabel',round(t(1:10:end)./sr))
                                    xlabel('Time, s')
                                end
                            end
                        end
                    end
                    if figs_flag
                        if save_wcc_fig == 1 % print to file
                            print(gcf,'-dpng','-r300','-loose',['wcc_' 'score' num2str(piecei) '_dim' num2str(traji) '_tr' num2str(triali) '_' datestr(now,'yymmdd-HHMMSS') '.png']);
                        else
                            pause % just inspect on the screen
                        end
                    end
            end
        end
        label_cc=[dataTrajs{traji},'_',method_flag];
        D{piecei}.(label_cc)=cor_vals(:,:,piecei);
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

%% CC Analysis - acceleration data
counter=0;
for piecei = 1:numel(D)
	% Preallocate vector to store correlation coefficients
    cor_vals=zeros(6*size(D{piecei}.A,2),1,numel(D));
        
    % Specify the same parameters that Andrew used
    counter=counter+1;
        
    switch method_flag
        case 'cc_and_gcorder'
            maxlag = morders(counter); %lag is the same as the model order for gc
            window=length(D{piecei}.A); %Whole length of the piece. Dobri suggests not to use this. Will discuss.
            overlap=0; %No overlap
        case 'wcc'
            sr = 100; % Hz. 
            maxlag=round(max_lag*sr); % 2 seconds
            window=round(win_len*sr); % 5 seconds 
            overlap=round(window/5); % half a window overlap
    end

	for triali=1:size(D{piecei}.A,2)
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
                            [c,l]=xcov(D{piecei}.A(row,:,triali),D{piecei}.A(col,:,triali),morders(counter),'coef'); %only thing I'm not sure about is 'coef' - normalizes the sequence 
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
                            [wcc,l,t]=corrgram(D{piecei}.A{triali}(row,:),D{piecei}.A{triali}(col,:),maxlag,window,overlap);
                            cor_vals(fcounter+6*(triali-1),1,piecei)=max(max(abs(wcc)));
                            if figs_flag
                                subplot(2,3,fcounter)
                                corrgram(D{piecei}.A{triali}(row,:),D{piecei}.A{triali}(col,:),maxlag,window,overlap)
                                xtickangle(30)
                                set(gca,'YTick',l(1:20:end))
                                set(gca,'YTickLabel',l(1:20:end)./sr)
                                ylabel('Lag, s')
                                set(gca,'XTick',t(1:8:end))
                                set(gca,'XTickLabel',round(t(1:8:end)./sr))
                                xlabel('Time, s')
                            end
                        end
                    end
                end
                if figs_flag
                    if save_wcc_fig == 1 % print to file
                        print(gcf,'-dpng','-r300','-loose',['wcc_' 'score' num2str(piecei) '_tr' num2str(triali) '_' datestr(now,'yymmdd-HHMMSS') '.png']);
                    else
                        pause % just inspect on the screen
                    end
                end
        end
    end
  	label_cc=['A_',method_flag];
   	D{piecei}.(label_cc)=cor_vals(:,:,piecei);
end
    
figure
for p=1:2
    subplot(2,2,1+(p-1)*2)
    boxplot(D{p}.A_wcc,reshape(meshgrid(1:8,1:6),1,[])')
end


%% Save data
if save_flag==1

    %Reconfigure cor_vals
    cor_vals_reconfig=[D{1}.X_processed_cc;D{2}.X_processed_cc];
    
    %Make vector for pair
    pair=repmat([1:6]',16,1);
    
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


    T=table(pair, cor_vals_reconfig,condition,trial,trial_collapsed,piece);
    filename='mada_cc3.xlsx';
    T.Properties.VariableNames = {'pair','cc', 'condition','trial', 'trial_collapsed','piece'};
    writetable(T,filename);
    
    plot(trial,cor_vals_reconfig)


end

