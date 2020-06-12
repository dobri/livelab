% load Desktop/DATA_star_2020_06.mat
% [WTCS,R,XWTS,W] = cwc_loop(lvtime,lvpos,censure);
% x [WTCS,R,XWTS,W] = cwc_loop(lvtime,lvposu,censure);
switch 4
    case 1
        CC = WTCS;
    case 2
        CC = R;
        %CC = mean(CC,4);
    case 3
        CC = XWTS;
        %CC = mean(XWTS,4);
        %CC = max(XWTS,[],4);
    case 4
        CC = log(abs(W));
        %CC = max(CC,[],4);
        CC = mean(CC,4);
        CC = CC-min(CC(:));
        CC = CC./max(CC(:));
end

close all


for lv=1:size(CC,4)
    figure(1)
    for tr=1:8
        subplot(4,4,tr+(lv-1)*8)
        imagesc(CC(:,:,tr,lv),[0 1])
    end
    colormap hot
end
pause(.2)


CClong = [];
for lv=1:size(CC,4)
    for tr=1:8
        cctemp = CC(:,:,tr,lv);
        cctemp = cctemp(~isnan(cctemp));
        cctemp(:,2) = tr;
        cctemp(:,3) = lv;
        CClong = vertcat(CClong,cctemp);
    end
end


for lv=1:size(CC,4)
    figure(2)
    subplot(1,2,lv)
    boxplot(CClong(CClong(:,3)==lv,1),CClong(CClong(:,3)==lv,2))
end
pause(.2)


cm = CM(DistMat>0);
dm = DistMat(DistMat>0);
for lv=1:size(CC,4)
    figure(3)
    for tr=1:8
        subplot(4,4,tr+(lv-1)*8)
        cc=CC(:,:,tr,lv);
        cc=cc(~isnan(cc));
        b=[dm*0+1 dm]\cc;
        [r,p] = corr(cc,dm);
        scatter(dm,cc)
        hold on
        plot(sort(dm),b(1)*(dm*0+1)+b(2)*sort(dm),'-or')
        hold off
        ylim([0 1])
        xlim([0 6.1])
        text(.1,.25,num2str([b',r,p,p<(.05/8)]','%10.4f'),'unit','normalized','fontsize',16)
    end
end
pause(.2)

for lv=1:size(CC,4)
    figure(4)
    for tr=1:8
        subplot(4,4,tr+(lv-1)*8)
        cc=CC(:,:,tr,lv);
        cc=cc(~isnan(cc));
        boxplot(cc,cm)
        [h,p,~,stats] = ttest2(cc(cm==0),cc(cm==1));
        disp([stats.tstat stats.df p])
        ylim([0 1])
    end
end
pause(.2)


for lv=1:size(CC,4)
    for tr=1:4
        cc1 = reshape(CC(:,:,tr  ,lv),[],1);
        cc1 = cc1(~isnan(cc1));
        cc2 = reshape(CC(:,:,tr+4,1),[],1);
        cc2 = cc2(~isnan(cc2));
        
        figure(5)
        subplot(2,4,tr + (lv-1)*4)
        [c,n] = hist(cc1);plot(n,c);hold on
        [c,n] = hist(cc2);plot(n,c);hold off
        [h,p,~,stats] = ttest(cc1,cc2);
        disp([stats.tstat stats.df p])
        
    end
    figure(6)
    subplot(1,2,lv)
    boxplot(CClong(CClong(:,3)==lv,1),CClong(CClong(:,3)==lv,2))
end
