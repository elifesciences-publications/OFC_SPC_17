
clear classify* missclas*
enssizes=205;
niter = 1000;

%% first classify the probe data on the basis of
% all "BD" trials and first 2 "AC" trials

max_A_trial = 2; max_B_trial = 3;

% calculate the likelihood of misclassifcation
window = 2;

neuron_index = nindex1;% (1:size(site,1))>0; %first in well behaved animals

[misclass_AoverC,xxxx] = bo_misclass_wrapper(neuron_index,psth_AB_trials,...
    psth_CD_trials, steps,niter, enssizes, max_A_trial, max_B_trial,window);

[misclass_AoverC_shuffle,xxxx_s] = bo_misclass_wrapper_shuffle(neuron_index,psth_AB_trials,...
    psth_CD_trials, steps,niter, enssizes, max_A_trial, max_B_trial,window);


neuron_index = nindex2; %then in reversed or null behavior @ probe

[misclass_CoverA] = bo_misclass_wrapper(neuron_index,psth_AB_trials,...
    psth_CD_trials, steps,niter, enssizes, max_A_trial, max_B_trial,window);

[misclass_CoverA_shuffle,xxxx_s] = bo_misclass_wrapper_shuffle(neuron_index,psth_AB_trials,...
    psth_CD_trials, steps,niter, enssizes, max_A_trial, max_B_trial,window);



% and then calculate a classification success rate relative to shuffled
% data

clear misclass_pct_above_shuffle*

for jj = size(misclass_CoverA,2):-1:1
    for ii = size(misclass_CoverA,3):-1:1
        
        misclass_pct_above_shuffle_g(jj,ii) = sum(misclass_AoverC_shuffle(:,jj,ii)>mean(misclass_AoverC(:,jj,ii),1));
        
        misclass_pct_above_shuffle_b(jj,ii) = sum(misclass_CoverA_shuffle(:,jj,ii)>mean(misclass_CoverA(:,jj,ii),1));
        
    end
end

%% do some plotting for the classification probe test data
% first intialize colors, index variables and y data to plot for good and
% bad animals relative to pvalue threshold

itr=1;
newstep = steps(2:end-1);
clr = [0 0 0;.5 .5 .5];

errorbins=130:10:200;

gsort = sort(squeeze((misclass_AoverC_shuffle(:,itr,:)))/max_B_trial);
bsort = sort(squeeze((misclass_CoverA_shuffle(:,itr,:)))/max_B_trial);

y_good = squeeze(nanmean(misclass_AoverC(:,itr,:),1)/max_B_trial);
ind_g = squeeze(misclass_pct_above_shuffle_g(itr,:)<.05*niter);
ind_g(sum(squeeze(misclass_AoverC_shuffle(:,itr,:))/max_B_trial>.5)<(niter*.05))=0;

y_bad = squeeze(nanmean(misclass_CoverA(:,itr,:),1)/max_B_trial);
ind_b = squeeze(misclass_pct_above_shuffle_b(itr,:)<.05*niter);
ind_b(sum(squeeze(misclass_CoverA_shuffle(:,itr,:))/max_B_trial>.5)<(niter*.05))=0;


xlabel('time from cue onset'),ylabel('classification accuracy (B=A; D=C)')


%% lineplot - patch bounds == 95%centile
figure;  subplot(1,1,1),  hold on;
errorbins=130:1:200;

mean_data = squeeze(nanmean(misclass_AoverC_shuffle(:,itr,errorbins),1)/max_B_trial)';
variance_data = squeeze(sort(misclass_AoverC_shuffle(:,itr,errorbins))/(max_B_trial));
numsessions_per_pct = floor(niter*.05);
upperbound = variance_data(end-numsessions_per_pct,:);
numsessions_per_pct = ceil(niter*.05);

lowerbound = variance_data(numsessions_per_pct,:);

h=patch([newstep(errorbins),newstep(errorbins(end:-1:1))],[lowerbound,upperbound(end:-1:1)],clr(1,:)+.3);hold on,
set(h,'LineStyle','none')
xlim([-5 10]),ylim([0 1]),set(gca,'TickDir','out','LineWidth',1),box off
plot(newstep,y_good,'linewidth',2.5,'color',clr(1,:)),hold on
scatter(newstep(ind_g),y_good(ind_g),150,clr(1,:),'o','filled'),hold on
plot(newstep,squeeze(nanmean(misclass_AoverC_shuffle(:,itr,:),1)/max_B_trial),':','color',clr(1,:))

errorbins=errorbins+1;
xlabel('time from cue onset'),ylabel('classification accuracy (B=A; D=C)')
subplot(1,1,1),hold on;


mean_data = squeeze(nanmean(misclass_CoverA_shuffle(:,itr,errorbins),1)/max_B_trial)';
variance_data = squeeze(sort(misclass_CoverA_shuffle(:,itr,errorbins))/(max_B_trial));
numsessions_per_pct = floor(niter*.05);
upperbound = variance_data(end-numsessions_per_pct,:);
numsessions_per_pct = ceil(niter*.05);

lowerbound = variance_data(numsessions_per_pct,:);

h=patch([newstep(errorbins),newstep(errorbins(end:-1:1))],[lowerbound,upperbound(end:-1:1)],clr(2,:)+.3);hold on,
set(h,'LineStyle','none')

xlim([-5 10]),ylim([0 1]),set(gca,'TickDir','out','LineWidth',1),box off

plot(newstep,y_bad,'linewidth',2.5,'color',clr(2,:)),hold on
scatter(newstep(ind_b),y_bad(ind_b),150,clr(2,:),'o','filled'),hold on
plot(newstep,squeeze(nanmean(misclass_CoverA_shuffle(:,itr,:),1)/max_B_trial),':','color',clr(2,:))
scatter(newstep(ind_g),y_good(ind_g),150,clr(1,:),'o','filled'),hold on
plot(newstep,squeeze(nanmean(misclass_AoverC_shuffle(:,itr,:),1)/max_B_trial),':','color',clr(1,:))
plot(newstep,y_good,'linewidth',2.5,'color',clr(1,:)),hold on


%% HISTOGRAMS
% plot the histogram of shuffled data relative to classification mean for a
% specific test bin for inset panels on figure 6

clear test*

testbin=188;
pltct=0;
pltct=pltct+2;

test_stat(:,1) = misclass_AoverC_shuffle(:,itr,testbin)/3;%-mean(misclass_AoverC(:,4,190));
test_stat(:,2) = misclass_CoverA_shuffle(:,itr,testbin)/3;%-mean(misclass_CoverA(:,4,190));

figure(12341);
histvec=0.0:.15:1.15;

subplot(3,2,pltct-1),h=bar(histvec+.06,histc(test_stat(:,1),histvec));
shading flat,box off,xlim([-.150 1.15]),ylim([0 500]),set(h,'barWidth',.8,'facecolor',[.2 .2 .2]+.1)
hold on,text(.75,300,['p = ',num2str((sum(test_stat(:,1)>mean(misclass_AoverC(:,itr,testbin)/3),1)+1)/(niter+1))])
hold on,text(0+.75,360,'A > C')
xlabel('accuracy (B=A; D=C)')


plot ([mean(misclass_AoverC(:,itr,testbin)/3) mean(misclass_AoverC(:,itr,testbin)/3)],[0 250],'r--')

subplot(3,2,pltct),h=bar(histvec+.06,histc(test_stat(:,2),histvec));
shading flat,box off,xlim([-.150 1.15]),ylim([0 500]),set(h,'barWidth',.8,'facecolor',[.2 .2 .2]+.1)
hold on,text(.75,300,['p = ',num2str((sum(test_stat(:,2)>mean(misclass_CoverA(:,itr,testbin)/3),1)+1)/(niter+1))])
hold on,text(0+.75,360,'C >/= A')
plot ([mean(misclass_CoverA(:,itr,testbin)/3) mean(misclass_CoverA(:,itr,testbin)/3)],[0 250],'r--')
xlabel('accuracy (B=A; D=C)')




