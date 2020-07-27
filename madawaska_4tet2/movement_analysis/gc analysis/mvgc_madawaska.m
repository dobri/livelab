%% mvgc_madawaska
% 
% This script completes a granger-causality analysis of the Madawaska data
% To use it, load in the variable 'D' from the MASTER_preprocess.m script 

% Before starting, make sure the MVGC toolbox is added to your path
addpath(genpath('~/Desktop/MATLAB/toolboxes/mvgc_v1.0'))

% Load data
%load('D.mat')

% Flag for saving data. Set to 1 if you want this loop to save a
% spreadsheet of the data. 0 if no

save_flag=1;

% Get fieldnames
dataTrajs=fieldnames(D{1})';

for piecei = 1:numel(D)
    for traji = 1:numel(fieldnames(D{piecei}))
        %% Parameters
        X=D{piecei}.(dataTrajs{traji});

        ntrials   = size(X,3);     % number of trials
        nobs      = size(X,2);   % number of observations per trial
        nvars     = size(X,1);      % number of variables
        
        bsize     = [];     % permutation test block size: empty for automatic (uses model order)

        regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
        icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

        morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
        momax     = 20;     % maximum model order for model order estimation

        acmaxlags = [];   % maximum autocovariance lags (empty for automatic calculation)

        tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
        alpha     = 0.05;   % significance level for significance test
        mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

        fs        = 7.6935;    % sample rate (Hz)
        fres      = [];     % frequency resolution (empty for automatic calculation)

        seed      = 0;      % random seed (0 for unseeded)

        %% Model order estimation 

        % Preallocate some vectors

        AIC_matrix=zeros(momax,ntrials);
        moAIC_matrix=zeros(1,ntrials);

        % Loop through all trials

        for triali=1:size(X,3)

        % Calculate information criteria up to specified maximum model order.

            [AIC,~,moAIC,~] = tsdata_to_infocrit(X(:,:,triali),momax,icregmode);
            AIC_matrix(:,triali)=AIC;

            moAIC_matrix(1,triali)=moAIC;

            % Plot information criteria.

            %figure(1); clf;
            %plot_tsdata([AIC]',{'AIC'},1/fs);
            %title('Model order estimation');
            %pause

        end

        morder=max(moAIC_matrix);
        label_morder=[dataTrajs{traji},'_morder'];
        D{piecei}.(label_morder)=morder;
        % This is the model order I'll use for every trial

        %% VAR model estimation (<mvgc_schema.html#3 |A2|>)

        % Preallocate some vectors
        GC_data=zeros(nvars,nvars,ntrials);
        pval_data=zeros(nvars,nvars,ntrials);
        sig_data=zeros(nvars,nvars,ntrials);
        cd_data=zeros(1,ntrials);

        % Loop through all trials

        for triali=1:size(X,3)

            % Estimate VAR model of selected order from data.

            [A,SIG] = tsdata_to_var(X(:,:,triali),morder,regmode);

            % Check for failed regression

            assert(~isbad(A),'VAR estimation failed');

            % Autocovariance calculation 
            % we calculate the autocovariance sequence G according to the
            % VAR model, to as many lags as it takes to decay to below the numerical
            % tolerance level, or to acmaxlags lags if specified (i.e. non-empty).

            [G,info] = var_to_autocov(A,SIG,acmaxlags);

            % Check for errors

            var_info(info,true); % report results (and bail out on error)

            % Calculate time-domain pairwise-conditional causalities 

            F = autocov_to_pwcgc(G);
            GC_data(:,:,triali)=F;

            % Check for failed GC calculation

            assert(~isbad(F,false),'GC calculation failed');

            % Significance test using theoretical null distribution, adjusting for multiple
            % hypotheses.

            pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
            pval_data(:,:,triali)=pval;
            label_pval=[dataTrajs{traji},'_pval'];
            D{piecei}.(label_pval)=pval_data;

            sig  = significance(pval,alpha,mhtc);
            sig_data(:,:,triali)=sig;

            % Plot time-domain causal graph, p-values and significance.

            figure(2); clf;
            subplot(1,3,1);
            plot_pw(F);
            title('Pairwise-conditional GC');
            subplot(1,3,2);
            plot_pw(pval);
            title('p-values');
            subplot(1,3,3);
            plot_pw(sig);
            title(['Significant at p = ' num2str(alpha)])

            % Calculate Seth's causal density (cd) measure, the mean pairwise-conditional causality. 

            cd = mean(F(~isnan(F)));
            cd_data(1,triali)=cd;
            label_gc=[dataTrajs{traji},'_gc'];
            D{piecei}.(label_gc)=GC_data;  
        end
        
    end
end

%% Save data
if save_flag == 1
                
	% Make table of the raw gc scores for each pair and save
	GCdata_reconfig1=[];
	GCdata_reconfig2=[];

	for piecei = 1:numel(D)
        for zz=1:size(D{piecei}.gc_data,3)
            for yy=1:size(D{piecei}.gc_data,2)
                for xx=1:size(D{piecei}.gc_data,1)
                    step=D{piecei}.gc_data(xx,yy,zz);
                    if piecei==1
                        GCdata_reconfig1=[GCdata_reconfig1,step];
                    else
                        GCdata_reconfig2=[GCdata_reconfig2,step];
                    end
                end
            end
        end
	end

    GCdata_reconfig1=GCdata_reconfig1(~isnan(GCdata_reconfig1))';
    GCdata_reconfig2=GCdata_reconfig2(~isnan(GCdata_reconfig2))';
    GCdata_reconfig_all=[GCdata_reconfig1;GCdata_reconfig2];

    %Make vector correponding to mechanical vs expressive trials
    condition1=[1+zeros(12,1);2+zeros(12,1);1+zeros(12,1);2+zeros(12,1);...
        1+zeros(12,1);2+zeros(12,1);1+zeros(12,1);2+zeros(12,1)];
    condition2=flip(condition1);
    condition=[condition1;condition2];

    %Make vector for time (trial)
    trial=[1+zeros(12,1);2+zeros(12,1);3+zeros(12,1);4+zeros(12,1);...
        5+zeros(12,1);6+zeros(12,1);7+zeros(12,1);8+zeros(12,1);];
    trial=repmat(trial,2,1);

    %Make vector for time collapsed across mechnical vs. expressive
    trial_collapsed = repmat(repelem([1:4]',24),2,1)

    %Make vector for piece
    piece=repelem([1;2],96);

    %Make vector for pair
    pair=repmat([1:12]',16,1);

    T=table(pair,GCdata_reconfig_all, condition, trial, trial_collapsed, piece);
    T.Properties.VariableNames = {'pair','gc', 'condition','trial', 'trial_collapsed','piece'}
    filename='mada_gc3.xlsx';
    writetable(T,filename);
end





