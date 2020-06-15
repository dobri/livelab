%{
addpath('~/Desktop/livelab/star/')
addpath('~/Desktop/livelab/star/wavelet-coherence/')
load Desktop/DATA_star_2020_06.mat
[WTCS,R,XWTS,W] = cwc_loop(lvtime,lvpos,censure);
save('~/Desktop/livelab/star/xwt_summary.mat','WTCS','R','XWTS','W','CC')
%}
% 
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
        % sum(sum(sum(mean(XWTS,4)==0))) ans = 0
    case 4
        CC = log(abs(W));
        %CC = log(abs(W));
        %CC = CC-min(CC(:));
        %CC = CC./max(CC(:));
        CC(:,:,:,3) = mean(CC,4);
end

close all


for lv=1:size(CC,4)
    figure(1)
    for tr=1:8
        subplot(6,4,tr+(lv-1)*8)
        imagesc(CC(:,:,tr,lv),[0 7])
        axis square
    end
    colormap hot
end
pause(.2)


CClong = [];
for lv=1:size(CC,4)
    for tr=1:8
        cctemp = CC(:,:,tr,lv);
        cctemp = cctemp(~isnan(cctemp));
        cctemp(:,2) = (1:numel(cctemp));
        cctemp(:,3) = tr;
        cctemp(:,4) = lv;
        CClong = vertcat(CClong,cctemp);
    end
end
CClong(:,5) = sum(CClong(:,3)==[5 6 7 8],2);
CClong(:,6) = sum(CClong(:,3)==[3 4 7 8],2);
CClong(:,7) = sum(CClong(:,3)==[2 4 6 8],2);
%{
fid=fopen(['xwt_long.csv'],'w');
fprintf(fid,'%s,','C','pp','tr','lv','eyes','groove','tempo');fprintf(fid,'\n');
for r=1:size(CClong,1);fprintf(fid,'%8.4f,',CClong(r,:));fprintf(fid,'\n');end
fclose(fid);
%}

for lv=1:size(CC,4)
    figure(2)
    subplot(1,3,lv)
    boxplot(CClong(CClong(:,4)==lv,1),CClong(CClong(:,4)==lv,3))
end
pause(.2)


if 0
    cm_type=4;
    cm = CM(:,:,cm_type);
    cm = cm(DistMat>0);
    dm = DistMat(DistMat>0);
    for lv=1:size(CC,4)
        figure(3)
        for tr=1:8
            subplot(6,4,tr+(lv-1)*8)
            cc=CC(:,:,tr,lv);
            cc=cc(~isnan(cc));
            b=[dm*0+1 dm]\cc;
            [r,p] = corr(cc,dm);
            scatter(dm,cc)
            hold on
            plot(sort(dm),b(1)*(dm*0+1)+b(2)*sort(dm),'-or')
            hold off
            %ylim([0 1])
            xlim([0 6.1])
            text(.1,.25,num2str([b',r,p,p<(.05/8)]','%10.4f'),'unit','normalized','fontsize',16)
        end
    end
    pause(.2)
    
    
    for lv=1:size(CC,4)
        figure(4)
        for tr=1:8
            subplot(6,4,tr+(lv-1)*8)
            cc=CC(:,:,tr,lv);
            cc=cc(~isnan(cc));
            boxplot(cc,cm)
            [h,p,~,stats] = ttest2(cc(cm==0),cc(cm==1));
            disp([stats.tstat stats.df p])
            %ylim([0 1])
        end
    end
    pause(.2)
end


for lv=1:size(CC,4)
    for tr=1:4
        cc1 = reshape(CC(:,:,tr  ,lv),[],1);
        cc1 = cc1(~isnan(cc1));
        cc2 = reshape(CC(:,:,tr+4,1),[],1);
        cc2 = cc2(~isnan(cc2));
        
        figure(5)
        subplot(3,4,tr + (lv-1)*4)
        [c,n] = hist(cc1);plot(n,c);hold on
        [c,n] = hist(cc2);plot(n,c);hold off
        [h,p,~,stats] = ttest(cc1,cc2);
        disp([stats.tstat stats.df p])
        
    end
    figure(6)
    subplot(1,3,lv)
    boxplot(CClong(CClong(:,4)==lv,1),CClong(CClong(:,4)==lv,3))
end


% CC to SPSS
Cmat_spss=[];
for lv=1:size(CC,4)
    Cmat=[];
    for c=1:1
        for trial=1:8
            temp = CC(:,:,trial,lv);
            temp(isnan(temp))=0;
            Cmat(:,trial)=mean(temp+temp')';
        end
        Cmat_spss=vertcat(Cmat_spss,horzcat(ones(size(Cmat,1),1)*lv,Cmat));
    end
    %{
    fid=fopen(['xwt_lv' num2str(lv,'%1.0f') '.csv'],'w');
    fprintf(fid,'%s,','000','001','010','011','100','101','110','111');fprintf(fid,'\n');
    for r=1:size(Cmat,1);fprintf(fid,'%8.4f,',Cmat(r,:));fprintf(fid,'\n');end
    fclose(fid);
    %}
end
