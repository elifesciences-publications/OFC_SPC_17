function [output,output2] = bo_ttest_cueresponse_inhib(psth_AB_background,psth_AB_trials,psth_CD_background,psth_CD_trials,indext,steps)
indext = logical(indext);
psth_AB_ttest = psth_AB_trials(indext,:,:,:);
psth_CD_ttest = psth_CD_trials(indext,:,:,:);
psth_AB_ttestback = psth_AB_background(indext,:,:,:);
psth_CD_ttestback = psth_CD_background(indext,:,:,:);


cueons=[0 0 2.5];
cueoffs=[10 2.5 10];

for ii=1
cueon = cueons(ii);
cueoff = cueoffs(ii);

cuetime = find(steps>cueon,1,'first'):find(steps<cueoff,1,'last');
for i = sum(indext):-1:1
    for j = 1:4
        if j<2.5
                [sign_p(i,j),sign_h(i,j)] = ranksum(nanmean(squeeze(psth_AB_ttest(i,j,1:6,cuetime)),2),squeeze(nanmean(psth_AB_ttestback(i,j,1:6,:),4)),'tail','left','alpha',0.05);
            

        else
                [sign_p(i,j),sign_h(i,j)] = ranksum(nanmean(squeeze(psth_CD_ttest(i,j,1:6,cuetime)),2),squeeze(nanmean(psth_CD_ttestback(i,j,1:6,:),4)),'tail','left','alpha',0.05);

        end
    end
end
%%

figure;hold on,
subplot(1,2,1);bar(sum(sign_h(:,:))/size(sign_h,1)),xlim([.5 4.5]),ylim([0 .5])
[xind,yind]=sort(-1.1*sign_h(:,3)-sign_h(:,4));
subplot(1,2,2);imagesc((sign_h(yind,:))),colormap('gray'),caxis([-.5 1.5])

%the fraction of neurons responding to one of the cues
sum(sum(sign_h(:,:),2)>0);

%%the fraction of neurons responding to each of the cues
output=(sum(sign_h(:,:),1))/size(sign_h,1);%/sum(sum(sign_h(:,:),1));
output2=sign_h;

howmany_ofAny_ABCD=sum(sum(sign_h(:,:),2)>0)%/sum(sum(sign_h(:,:),1));
howmany_ofAandB = sum(sum(sign_h(:,[1 3]),2)>1)
howmany_ofAorB = sum(sum(sign_h(:,[1 3]),2)>0)
howmany_ofCandD = sum(sum(sign_h(:,[2 4]),2)>1)
howmany_ofCorD = sum(sum(sign_h(:,[2 4]),2)>0)
howmany_ofAandD = sum(sum(sign_h(:,[1 4]),2)>1)
howmany_ofAorD = sum(sum(sign_h(:,[1 4]),2)>0)
howmany_ofCandB = sum(sum(sign_h(:,[2 3]),2)>1)
howmany_ofCorB = sum(sum(sign_h(:,[2 3]),2)>0)
howmany_withinpairsig = sum(((sum(sign_h(:,[1 3]),2)>1)+(sum(sign_h(:,[2 4]),2)>1))>0)
howmany_withoutpairsig = sum(((sum(sign_h(:,[2 3]),2)>1)+(sum(sign_h(:,[1 4]),2)>1))>0)

figure;hold on,
subplot(1,2,1);bar(sum(sign_h(:,:))/size(sign_h,1)),xlim([.5 4.5]),ylim([0 .5])
[xind,yind]=sort(-1.1*sign_h(:,3)-sign_h(:,4));
subplot(1,2,2);imagesc((sign_h(yind,:))),colormap('gray'),caxis([-.5 1.5])
end 

