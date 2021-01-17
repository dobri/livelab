%% CC_madawaska

% First, load in the data variable 'D'

% Load data
%load('D.mat')

% Flag for saving data. Set to 1 if you want this loop to save a
% spreadsheet of the data. 0 if no

save_flag=0;
figs_flag=1;
save_wcc_fig=0;
if save_wcc_fig==1
    close all
    figure('Position',[100 100 800*6 800])
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
    sr8 = 8; % Hz. After downsampling?
    sr100 = 100;
    win_len = 5; % seconds
    max_lag = 2; % seconds
end
num_trials = numel(D{1}.A);

if strcmp(method_flag,'cc_and_gcorder')
    morders=[D{1}.X_processed_morder,D{1}.X_detrended_processed_morder,...
    D{2}.X_processed_morder,D{2}.X_detrended_processed_morder];
end

dataTrajs={'X_processed','X_detrended_processed','A'};

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
        
        % Preallocate vector to store correlation coefficients
        % cor_vals=zeros(6*size(D{piecei}.(dataTrajs{traji}),3),1,numel(D)); % there's something weird here.
        cor_vals = [];
        
        % Specify the same parameters that Andrew used
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
            switch 0
                case 0
                    %salient_points = (30:30:max(t/sr))';
                    salient_points = cumsum(randi(10,10,1)+30-5);
                case 1
                    % Insert here a line to read a clicktrack from a text file such as the marked labels exported from Audacity.
                    % salient_points = csvread....
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
                                if traji==3
                                    x = D{piecei}.A{triali}(row,:)';
                                    y = D{piecei}.A{triali}(col,:)';
                                else
                                    x = D{piecei}.(dataTrajs{traji})(row,:,triali);
                                    y = D{piecei}.(dataTrajs{traji})(col,:,triali);
                                end
                                [wcc,l,t]=corrgram(x,y,maxlag,window,overlap);
                                cor_vals(fcounter+6*(triali-1),1,piecei)=max(max(abs(wcc)));
                                if figs_flag
                                    % This makes a pseudo-time vector that 
                                    % runs in units of time defined by the 
                                    % time points in the salient_points vector.
                                    tn = time_interpolate(t./sr,salient_points); 

                                    %corrgram(D{piecei}.(dataTrajs{traji})(row,:,triali),D{piecei}.(dataTrajs{traji})(col,:,triali),maxlag,window,overlap)
                                    subplot(num_trials,6,fcounter+6*(triali-1))
                                    imagesc(tn,l,flip(wcc),[-1 1]);colorbar
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
                                    set(gca,'XTick',1:numel(salient_points))
                                    %set(gca,'XTick',t(1:10:end))
                                    %set(gca,'XTickLabel',round(t(1:10:end)./sr))
                                    %if triali==size(D{piecei}.(dataTrajs{traji}),3)
                                    %    xlabel('Time [seconds or sections]')
                                    %end
                                    if triali==1
                                        title(['Pair ' num2str(row) '-' num2str(col)])
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
        D{piecei}.(label_cc)=cor_vals(:,:,piecei);
    end
end

if save_wcc_fig == 1
    set(gcf,'visible','on')
    close all
end

return 

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
                        print(gcf,'-dpng','-r100','-loose',['wcc_' 'score' num2str(piecei) '_tr' num2str(triali) '_' datestr(now,'yymmdd-HHMMSS') '.png']);
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

