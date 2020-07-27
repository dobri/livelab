function M = create_matrix_for_mvgc(DATA, dataTraj, plotting_flag)

% This function takes the processed madawaska body sway data and gets it 
% in the correct form for the mvgc toolbox

% The input is the all the processed ensemble madawaska data from the 
% prepare_data_for_mvgc function.

% If you want to plot everything, set plotting_flag=1. If not, set it to 0.

% first, let's see how long each trial is
trialLengths=zeros(1,length(DATA));
for i=1:length(DATA)
    trialLengths(1,i)=length(DATA{i}.(dataTraj));
end

minVal=min(trialLengths);
% ok I'm going to cap all of them at minVal (2513) because that's the shortest

% now create the matrix M for MVGC
nvars=size(DATA{1}.(dataTraj),3);
nobs=minVal;
ntrials=length(DATA);

M=zeros(nvars,nobs,ntrials);

for tri=1:ntrials
    AP=permute(DATA{tri}.(dataTraj),[3,1,2]);
    M(:,:,tri)=AP(:,1:minVal);
end

%% Optional - plot the matrix to make sure it looks ok 
if plotting_flag ==1

    for i=1:size(M,3) %loop through all 8 trials
        for j=1:size(M,1) % loop through all 4 musicians
            plot(M(j,:,i))
            disp([i, j])
            pause
        end
    end
end

