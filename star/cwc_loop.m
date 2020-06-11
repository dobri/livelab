function CC = cwc_loop(t,data,censure)

ds=4;

plotting = 2;
saving = 1;

rows = meshgrid(1:33,1:33)';
cols = meshgrid(1:33,1:33);

source_target_pairs1 = rows(:);
source_target_pairs2 = cols(:);

CC = nan(33,33,8,2);
total_runs = 8*33*33*2;

for tr = 1:8
    data{tr}(:,:,censure)=[];
end

[b,a] = butter(2,.2/50,'high');
for tr=1:8
    for lv=1:2
        for pp=1:33
            % Just in case clean a little
            data{tr}(data{tr}(:,lv,pp)==-inf,lv,pp)=0;
            data{tr}(data{tr}(:,lv,pp)==inf,lv,pp)=0;
            data{tr}(isnan(data{tr}(:,lv,pp)),lv,pp)=0;
            data{tr}(:,lv,pp) = filtfilt(b,a,data{tr}(:,lv,pp));
            %data{trial}(:,lv,pp)=smooth(data{trial}(:,lv,pp),10);
        end
    end
end


runs_n=0;
tic
for lv=1:2
    for tr=1:8
        CCtemp=nan(33*33,1);
        %for st=1:33^2
        parfor st=1:33^2
            row=source_target_pairs1(st);
            col=source_target_pairs2(st);
            if row>col
                if saving == 1
                    pic_name = ['~/cwc_pics/cwc_pps-' num2str(row,'%02.f') '-' num2str(col,'%02.f') '_tr' num2str(tr,'%02.f') '_lv' num2str(lv,'%02.f')];
                else
                    pic_name = [];
                end
                CCtemp(st) = cwc_test([t{tr}(1:ds:end) data{tr}(1:ds:end,lv,row)],...
                    [t{tr}(1:ds:end) data{tr}(1:ds:end,lv,col)],...
                    plotting,pic_name);
            end
        end
        for st=1:33^2
            runs_n=runs_n+1;
            row=source_target_pairs1(st);
            col=source_target_pairs2(st);
            CC(row,col,tr,lv) = CCtemp(st);
        end
        fprintf('That''s %.3f proportion of all runs completed in %.2f minutes.\n',runs_n/total_runs,toc/60)
    end
end


end