%% tackles all steps of correlation display
% however assumes that data is coming over from within bo_runmain


cues = 1:4;
clear rho*

rhotimeStart = 0; % start analysis window
rhotimeStop = 10; % stop analysis window


clr = [0 0 0;170 170 170;0,101,245;204,20,0]/255;
 figure;

clr =  [0.3 0.77 0.99;...
        0.88 0.88 0.14;...
        0.55 0.55 0.15;...
        0.4 0.4 0.8];...

meanpsth2=neural_resp(nindex,:,:);

clear cue_*
for nn = 1:size(meanpsth2,1)
    cue_mean(nn,cues)=squeeze(mean(meanpsth2(nn,:,find(steps>rhotimeStart,1,'first'):find(steps<rhotimeStop,1,'last')),3))/bin;
    cue_meanBK(nn,cues)=squeeze(mean(meanpsth2(nn,:,find(steps>-30,1,'first'):find(steps<-20,1,'last')),3))/bin;
    
end

sjsj = 1:size(meanpsth2,1);


subplot(4,2,1)
scatter(cue_mean(sjsj,1)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,3)-mean(cue_meanBK(sjsj,:),2),[],clr(1,:))

set(gca,'TickDir','out','LineWidth',1),box off
xlabel('hz to B over baseline'),ylabel('hz to A over baseline')

subplot(4,2,3)
scatter(cue_mean(sjsj,1)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,4)-mean(cue_meanBK(sjsj,:),2),[],clr(3,:))

set(gca,'TickDir','out','LineWidth',1),box off
xlabel('hz to B over baseline'),ylabel('hz to C over baseline')

subplot(4,2,2)
scatter(cue_mean(sjsj,2)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,3)-mean(cue_meanBK(sjsj,:),2),[],clr(2,:))

set(gca,'TickDir','out','LineWidth',1),box off
xlabel('hz to D over baseline'),ylabel('hz to A over baseline')

subplot(4,2,4)
scatter(cue_mean(sjsj,2)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,4)-mean(cue_meanBK(sjsj,:),2),[],clr(4,:))

set(gca,'TickDir','out','LineWidth',1),box off
xlabel('hz to D over baseline'),ylabel('hz to C over baseline')

ctype = 'pearson';
[rho(1), rhop(1)]=corr(cue_mean(sjsj,1)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,3)-mean(cue_meanBK(sjsj,:),2),'type',ctype);
[rho(2), rhop(2)]=corr(cue_mean(sjsj,1)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,4)-mean(cue_meanBK(sjsj,:),2),'type',ctype);
[rho(3), rhop(3)]=corr(cue_mean(sjsj,2)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,3)-mean(cue_meanBK(sjsj,:),2),'type',ctype);
[rho(4), rhop(4)]=corr(cue_mean(sjsj,2)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,4)-mean(cue_meanBK(sjsj,:),2),'type',ctype);


for ii = 1:4
subplot(4,2,ii)
hold on,
text(-1,28.2,['r = ',num2str(rho(ii))]),
text(-1,21.6,['p = ',num2str(rhop(ii))]);
xlim([-15 35.5]),ylim([-20 35.5])
end

%% then get figure 5D/5E histogram of one bootstrap; 


jjct=0;


for jj = 5:5:size(meanpsth2,1),
    jjct = jjct+1;
for ii = 1:10
clear cue_mean
clear cue_meanBK

sjsj = randperm(size(meanpsth2,1),jj);
for nn = 1:size(meanpsth2,1)
    cue_mean(nn,cues)=squeeze(mean(meanpsth2(nn,:,find(steps==rhotimeStart,1,'first'):find(steps==rhotimeStop,1,'last')),3));
    cue_meanBK(nn,cues)=squeeze(mean(meanpsth2(nn,:,find(steps==-30,1,'first'):find(steps==-0,1,'last')),3));
    
end

[rho(1), rhop(1)]=corr(cue_mean(sjsj,1)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,3)-mean(cue_meanBK(sjsj,:),2),'type','pearson');
[rho(2), rhop(2)]=corr(cue_mean(sjsj,1)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,4)-mean(cue_meanBK(sjsj,:),2),'type','pearson');
[rho(3), rhop(3)]=corr(cue_mean(sjsj,2)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,3)-mean(cue_meanBK(sjsj,:),2),'type','pearson');
[rho(4), rhop(4)]=corr(cue_mean(sjsj,2)-mean(cue_meanBK(sjsj,:),2),cue_mean(sjsj,4)-mean(cue_meanBK(sjsj,:),2),'type','pearson');



rho_agg(jjct,ii,:) = rho;
end
end



subplot(4,2,[1:2]+4)
xhistvalues =-0.5:0.075:1;
yhistvalues=histc(squeeze(rho_agg(ceil(2*jjct/4),:,:)),xhistvalues);
bar(xhistvalues,yhistvalues,1.25,'edgecolor','none'),colormap(clr),xlim([-0.25,0.85])
set(gca,'TickDir','out','LineWidth',1),box off
jjvector=5:5:size(meanpsth2,1);
subplot(4,2,[3:4]+4)


for ii = 1:4
    errorbar(jjvector,squeeze(nanmean(rho_agg(:,:,ii),2))',2*squeeze(nanstd(rho_agg(:,:,ii),[],2))'/sqrt(sum(nindex>0)),'color',clr(ii,:),'linewidth',2),hold on
    xlabel('ensemble size'),ylabel('correlation between cues')
    set(gca,'TickDir','out','LineWidth',1),box off
    hold on
    plot([jjvector(ceil(2*jjct/4)),jjvector(ceil(2*jjct/4))],[-.5 1],'k--')
    ylim([-0.25,0.85])
    xlim([0 max(jjvector)])
end
