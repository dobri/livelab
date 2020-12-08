%% CC_madawaska

% First, load in the data variable 'D'

% Load data
%load('D.mat')

% Flag for saving data. Set to 1 if you want this loop to save a
% spreadsheet of the data. 0 if no

save_flag=0;

%Specify parameters
%Which method are we going to use? Dobri =0 or Andrew = 1?
method_flag=1;

%if method_flag
%end
%haven't tried these different methods yet. Keep doing Andrew's method for now

dataTrajs={'X_processed','X_detrended_processed'};

morders=[D{1}.X_processed_morder,D{1}.X_detrended_processed_morder,...
    D{2}.X_processed_morder,D{2}.X_detrended_processed_morder];


%% CC analysis
counter=[];

for piecei = 1:numel(D)
    
    for traji = 1:numel(dataTrajs)

        % Preallocate vector to store correlation coefficients
        cor_vals=zeros(6*size(D{piecei}.(dataTrajs{traji}),3),1,numel(D));
        
        % Specify the same parameters that Andrew used
        counter=counter+1;
        maxlag=morders(counter); %lag is the same as the model order for gc
        window=length(D{piecei}.(dataTrajs{traji})); %Whole length of the piece. Dobri suggests not to use this. Will discuss.
        overlap=0; %No overlap

        for triali=1:size(D{piecei}.(dataTrajs{traji}),3)
            
            %take the maximum unsigned CC coefficient for each of the 6 possible
            %pairs of musicians for each trial
            cor_vals(1+6*(triali-1),1,piecei)=max(abs(corrgram(D{piecei}.(dataTrajs{traji})(1,:,triali),D{piecei}.(dataTrajs{traji})(2,:,triali),maxlag, window,overlap)),[],'all');
            cor_vals(2+6*(triali-1),1,piecei)=max(abs(corrgram(D{piecei}.(dataTrajs{traji})(1,:,triali),D{piecei}.(dataTrajs{traji})(3,:,triali),maxlag,window,overlap)),[],'all');
            cor_vals(3+6*(triali-1),1,piecei)=max(abs(corrgram(D{piecei}.(dataTrajs{traji})(1,:,triali),D{piecei}.(dataTrajs{traji})(4,:,triali),maxlag,window,overlap)),[],'all');
            cor_vals(4+6*(triali-1),1,piecei)=max(abs(corrgram(D{piecei}.(dataTrajs{traji})(2,:,triali),D{piecei}.(dataTrajs{traji})(3,:,triali),maxlag,window,overlap)),[],'all');
            cor_vals(5+6*(triali-1),1,piecei)=max(abs(corrgram(D{piecei}.(dataTrajs{traji})(2,:,triali),D{piecei}.(dataTrajs{traji})(4,:,triali),maxlag,window,overlap)),[],'all');
            cor_vals(6+6*(triali-1),1,piecei)=max(abs(corrgram(D{piecei}.(dataTrajs{traji})(3,:,triali),D{piecei}.(dataTrajs{traji})(4,:,triali),maxlag,window,overlap)),[],'all');
            
        end
        
        label_cc=[dataTrajs{traji},'_cc'];
        D{piecei}.(label_cc)=cor_vals(:,:,piecei); 
        
    end
    
    
end

%save('D','D')


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

