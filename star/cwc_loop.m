function [CC1,RR,CC2,WW] = cwc_loop(t,data,censure)

ds=4;

plotting = 2;
saving = 1;

rows = meshgrid(1:33,1:33)';
cols = meshgrid(1:33,1:33);

source_target_pairs1 = rows(:);
source_target_pairs2 = cols(:);

CC1 = nan(33,33,8,2);
RR = nan(33,33,8,2);
CC2 = nan(33,33,8,2);
WW = nan(33,33,8,2);
total_runs = 8*33*33*2;

for tr = 1:8
    data{tr}(:,:,censure)=[];
end

[b,a] = butter(2,.3/50,'high');
for tr=1:8
    for lv=1:2
        for pp=1:33
            % Just in case clean a little
            data{tr}(data{tr}(:,lv,pp)==-inf,lv,pp)=0;
            data{tr}(data{tr}(:,lv,pp)==inf,lv,pp)=0;
            data{tr}(isnan(data{tr}(:,lv,pp)),lv,pp)=0;
            %data{tr}(:,lv,pp) = filtfilt(b,a,data{tr}(:,lv,pp));
            %data{trial}(:,lv,pp)=smooth(data{trial}(:,lv,pp),10);
        end
    end
end


runs_n=0;
tic
for lv=1:2
    for tr=1:8
        CC1temp=nan(33*33,1);
        CC2temp=nan(33*33,1);
        Rtemp=nan(33*33,1);
        Wtemp=nan(33*33,1);
        %for st=1:33^2
        parfor st=1:33^2
            row=source_target_pairs1(st);
            col=source_target_pairs2(st);
            if row>col
                if saving == 1
                    pic_name2 = ['~/cwc_pics/wtc_pps-' num2str(row,'%02.f') '-' num2str(col,'%02.f') '_tr' num2str(tr,'%02.f') '_lv' num2str(lv,'%02.f')];
                    pic_name3 = ['~/cwc_pics/xwt_pps-' num2str(row,'%02.f') '-' num2str(col,'%02.f') '_tr' num2str(tr,'%02.f') '_lv' num2str(lv,'%02.f')];
                else
                    pic_name2 = [];
                    pic_name3 = [];
                end
                [CC1temp(st),Rtemp(st),CC2temp(st),Wtemp(st)] = cwc_test([t{tr}(1:ds:end) data{tr}(1:ds:end,lv,row)],...
                    [t{tr}(1:ds:end) data{tr}(1:ds:end,lv,col)],...
                    plotting,pic_name2,pic_name3);
            end
        end
        for st=1:33^2
            runs_n=runs_n+1;
            row=source_target_pairs1(st);
            col=source_target_pairs2(st);
            CC1(row,col,tr,lv) = CC1temp(st);
            RR(row,col,tr,lv) = Rtemp(st);
            CC2(row,col,tr,lv) = CC2temp(st);
            WW(row,col,tr,lv) = Wtemp(st);
        end
        fprintf('That''s %.3f proportion of all runs completed in %.2f minutes.\n',runs_n/total_runs,toc/60)
    end
end


end