function time_out = time_interpolate(time_in,checkpoints)

checkpoints = checkpoints(checkpoints<time_in(end));

if size(time_in,2) > size(time_in,1)
    time_in = time_in';
end
time_out = time_in*0;
[~,ind] = min(abs(time_in-checkpoints(1)));
time_out(1:ind) = linspace(0,1,ind);
for tt=1:numel(checkpoints)-1
    [~,ind1] = min(abs(time_in-checkpoints(tt)));
    [~,ind2] = min(abs(time_in-checkpoints(tt+1)));
    time_out(ind1+1:ind2) = linspace(tt+1/(ind2-ind1),tt+1,ind2-ind1);
end
[~,ind] = min(abs(time_in-checkpoints(end)));
time_out(ind+1:end) = linspace(tt+1+1/(numel(time_in)-ind),tt+1+1,numel(time_in)-ind);