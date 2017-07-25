%% for ofc data piped from bo_runmain


clr = [0 0 0;170 170 170;0,101,245;204,20,0]/255;
     psthplot = (sgolayfilt(squeeze(nanmean((neural_resp(nindex,ii,:)),1)),1,3)');


figure;
endind=find(steps>0,1,'first');
endind2=find(steps<10,1,'last');

[maxval]=(nanmean(neural_resp(:,3,endind:endind2)-neural_resp(:,4,endind:endind2),3))./(.1+nanmean(neural_resp(:,3,endind:endind2)+neural_resp(:,4,endind:endind2),3));
[xval,xind]=sort(maxval);

subplot(1,4,1)
imagesc(steps(1:length(psthplot)),1:sum(nindex),squeeze(neural_resp(xind,3,1:length(psthplot))))
xlim([-10 25]),xlabel('time from cue B onset'),set(gca,'TickDir','out','LineWidth',1),box off
colormap('gray'),caxis([.5 .85])

subplot(1,4,2)
imagesc(steps(1:length(psthplot)),1:sum(nindex),squeeze(neural_resp(xind,4,1:length(psthplot))))
xlim([-10 25]),xlabel('time from cue D onset'),set(gca,'TickDir','out','LineWidth',1),box off
colormap('gray'),caxis([.5 .85])


subplot(1,4,3)
imagesc(steps(1:length(psthplot)),1:sum(nindex),squeeze(neural_resp(xind,1,1:length(psthplot))))
xlim([-10 25]),xlabel('time from cue A onset'),set(gca,'TickDir','out','LineWidth',1),box off
colormap('gray'),caxis([.5 1])


subplot(1,4,4)
imagesc(steps(1:length(psthplot)),1:sum(nindex),squeeze(neural_resp(xind,2,1:length(psthplot))))
xlim([-10 25]),xlabel('time from cue C onset'),set(gca,'TickDir','out','LineWidth',1),box off
colormap('gray'),caxis([.5 .85])


steps2=steps
nindex = logical(nindex)
nindex4 = xind(1:ceil(sum(nindex)/20));
nindex3 = xind(end-ceil(sum(nindex)/20):end);


neural_resp = AUC;
neural_resp=(medfilt1(neural_resp(nindex,:,:),3,[],3));


figure;
for ii = 1:4;
    hold on
    psthsub =  0;
    psthplot = ((squeeze(nanmean((neural_resp(nindex3,ii,:)),1)))');
    semplot = 2*((squeeze(nanstd((neural_resp(nindex3,ii,:)),1)))')/sqrt(length(nindex3));
    plot(steps2(1:length(psthplot)),psthplot-psthsub,'color',clr(ii,:),'lineWidth',2);
    patch([steps2(1:length(psthplot)) steps2(length(psthplot):-1:1)],[psthplot(1:end)+semplot psthplot(end:-1:1)-semplot]-psthsub,clr(ii,:),'linestyle','none')
    xlim([-15 20]),%ylim([.4 .8])
    xlabel('time from cue onset')
    ylabel('mean Hz')
    set(gca,'TickDir','out','LineWidth',1),box off
end

figure;

for ii = 1:4;
    hold on
    psthsub = 0;
    psthplot = ((squeeze(nanmean((neural_resp(nindex4,ii,:)),1)))');
    semplot = 2*((squeeze(nanstd((neural_resp(nindex4,ii,:)),1)))')/sqrt(length(nindex4));
    plot(steps2(1:length(psthplot)),psthplot-psthsub,'color',clr(ii,:),'lineWidth',2);
    patch([steps2(1:length(psthplot)) steps2(length(psthplot):-1:1)],[psthplot(1:end)+semplot psthplot(end:-1:1)-semplot]-psthsub,clr(ii,:),'linestyle','none')
    xlim([-15 20]),%ylim([.4 .8])
    xlabel('time from cue onset')
    ylabel('mean Hz')
    set(gca,'TickDir','out','LineWidth',1),box off
end
