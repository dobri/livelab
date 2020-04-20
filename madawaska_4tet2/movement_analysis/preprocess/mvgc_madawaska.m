%% mvgc_madawaska
% 
% This script completes a granger-causality analysis of the Madawaska data
% To use it, load in the matrix M from the script 'prepare_data_for_mvgc.m'

%% Parameters

ntrials   = 8;     % number of trials
nobs      = 2513;   % number of observations per trial
nvars     = 4;      % number of variables

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 7.6935;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)

X=M;

%% Model order estimation 

% Preallocate some vectors

AIC_matrix=zeros(momax,ntrials);
moAIC_matrix=zeros(1,ntrials);

% Loop through all trials

for triali=1:size(X,3)

% Calculate information criteria up to specified maximum model order.
    
    [AIC,~,moAIC,~] = tsdata_to_infocrit(M(:,:,triali),momax,icregmode);
    AIC_matrix(:,triali)=AIC;
    
    moAIC_matrix(1,triali)=moAIC;

    % Plot information criteria.

    %figure(1); clf;
    %plot_tsdata([AIC BIC]',{'AIC','BIC'},1/fs);
    %title('Model order estimation');
    %pause
    
end

morder=max(moAIC_matrix);
% This is the model order I'll use for every trial

%% VAR model estimation (<mvgc_schema.html#3 |A2|>)

% Preallocate some vectors

%A_matrix=zeros(momax,ntrials);
%SIG_matrix=zeros(1,ntrials);
%G_matrix
%info_matrix

GC_data=zeros(nvars,nvars,ntrials);
pval_data=zeros(nvars,nvars,ntrials);
sig_data=zeros(nvars,nvars,ntrials);
cd_data=zeros(1,ntrials);

% Loop through all trials

for triali=1:size(X,3)
    
    % Estimate VAR model of selected order from data.

    [A,SIG] = tsdata_to_var(M(:,:,triali),morder,regmode);

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

end


odd_trials=GC_data(:,:,[1,3,5,7]);
odd_pvals=pval_data(:,:,[1,3,5,7]);
odd_sig=sig_data(:,:,[1,3,5,7]);
odd_cd=cd_data(1,[1,3,5,7]);


    figure(1); clf;
    subplot(1,3,1);
    plot_pw(mean(odd_trials, 3));
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(mean(odd_pvals,3));
    title('p-values');
    subplot(1,3,3);
    plot_pw(mean(odd_sig,3));
    title(['Significant at p = ' num2str(alpha)])


even_trials=GC_data(:,:,[2,4,6,8]);
even_pvals=pval_data(:,:,[2,4,6,8]);
even_sig=sig_data(:,:,[2,4,6,8]);
even_cd=cd_data(1,[2,4,6,8]);


    figure(2); clf;
    subplot(1,3,1);
    plot_pw(mean(even_trials, 3));
    title('Pairwise-conditional GC');
    subplot(1,3,2);
    plot_pw(mean(even_pvals,3));
    title('p-values');
    subplot(1,3,3);
    plot_pw(mean(even_sig,3));
    title(['Significant at p = ' num2str(alpha)])


%% Save data

% Make a vector to correspond to odd or even trials

condition = [1;2;1;2;1;2;1;2];
trial=[1;2;3;4;5;6;7;8];

% Make table

T=table(cd_data',trial, condition);

% Save data
filename='madawaska.xlsx';
writetable(T,filename);




