%% CC_madawaska

% First, load in the matrix 'M'

%% Let's try Andrew's method first

% Preallocate a vector
cor_vals=zeros(6,size(M,3));
maxlag=17;
window=length(M);
overlap=0;

for triali=1:size(M,3)
    cor_vals(1,triali)=max(abs(corrgram(M(1,:,triali),M(2,:,triali),maxlag, window,overlap)),[],'all');
    cor_vals(2,triali)=max(abs(corrgram(M(1,:,triali),M(3,:,triali),maxlag,window,overlap)),[],'all');
    cor_vals(3,triali)=max(abs(corrgram(M(1,:,triali),M(4,:,triali),maxlag,window,overlap)),[],'all');
    cor_vals(4,triali)=max(abs(corrgram(M(2,:,triali),M(3,:,triali),maxlag,window,overlap)),[],'all');
    cor_vals(5,triali)=max(abs(corrgram(M(2,:,triali),M(4,:,triali),maxlag,window,overlap)),[],'all');
    cor_vals(6,triali)=max(abs(corrgram(M(3,:,triali),M(4,:,triali),maxlag,window,overlap)),[],'all');
end


mean(cor_vals)

cor_vals_reconfig=[cor_vals(:,1); cor_vals(:,2); cor_vals(:,3); cor_vals(:,4);...
    cor_vals(:,5);cor_vals(:,6);cor_vals(:,7);cor_vals(:,8);];
condition=[1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2]';
trial=[1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8]';

T=table(cor_vals_reconfig,condition,trial);
filename='mada_cc.xlsx';
writetable(T,filename);


%% Now let's try Dobri's suggestion
%sr=8;
%[c,l,t]=corrgram(M(1,:,1),M(2,:,1),round(sr*1),round(sr*5),round(sr*5*.5));


