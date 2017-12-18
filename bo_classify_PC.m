    load('probe_beh_vars.mat')

clear filename
filename(1) = {'PC1_OFC_collected_aligned_.mat'};
filename(2) = {'PC2_OFC_collected_aligned_.mat'};
bin = .25;

sortthreshold = 3; % remove neurons w/ deviation less than 3x noise

clr=[0 0 0;... %events 1 (black)
    .5 .5 .5;...  %events 2 (gray)
    .1 .2 1;...  %events 3 (blue)
    1 .1 .2];    %events 4 (red)

classifytype = 'trial_blocks'; %4 trials at a time, 6x
%classifytype = 'trial_by_trial'; %all 24 trials individually

ccount=0;

clear rho p
scatterindex=0;
for fi =1:2;
    
    %% collect and load the binned spike data to use in average cue response


    [site, neuronEnsNum, pokein, pokeout] = load_spc(filename{fi});
        
    clear sortstat
    for ii = size(site,1):-1:1, sortstat(ii)=site(ii,1).stats.sig2noise;end
    site = site(sortstat>sortthreshold,:);
    neuronEnsNum = neuronEnsNum(sortstat>sortthreshold,:);
    
    goodEns =  ismember(neuronEnsNum,beh_ens(beh_order(beh_num>0)));

    
    %site=site(goodEns,:);
    %neuronEnsNum = neuronEnsNum(goodEns,:);
    
    nindex = 1:size(site,1);
    
    [meanpsth, bsubpsth, AUC, steps, psth_AB_trials,psth_CD_trials,psth_AB_background,psth_CD_background]= bo_psthify_pc(site,neuronEnsNum,bin);
    
    cues = 1:4;
    
    
    %% extract the mean response to each cue
    cue_mean = zeros(size(site,1),4);   cue_meanBK = zeros(size(site,1),4);
    bsubCuemean = zeros(size(site,1),4);
    cue_mean_tr = zeros(size(site,1),24); cue_meanBK_tr = zeros(size(site,1),24);
    
    for trials=6:-1:1;
        for i = size(site,1):-1:1
            l = 0;
            for j = cues
                if j<2.5
                    meanpsth(i,j,:) =  squeeze(nanmean(psth_AB_trials(i,j,trials,:),3));
                    
                else
                    meanpsth(i,j,:) =  squeeze(nanmean(psth_CD_trials(i,j,trials,:),3));
                    
                end
            end
        end
        
        
        for nn = 1:size(site,1)
            
            cue_mean_tr(nn,4*(trials-1)+1:4*(trials)) = squeeze(nanmean(meanpsth(nn,:,find(steps>0,1,'first'):find(steps<10,1,'last')),3));%-squeeze(mean(mean(yyy(nn,:,:,:),4),3));
            cue_meanBK_tr(nn,4*(trials-1)+1:4*(trials)) = squeeze(nanmean(meanpsth(nn,:,find(steps>-30,1,'first'):find(steps<-20,1,'last')),3));%-squeeze(mean(mean(yyy(nn,:,:,:),4),3));
            
        end
    end
    
    cue_mean_tr(isnan(cue_mean_tr))=0;
    cue_meanBK_tr(isnan(cue_mean_tr))=0;
    %% get ready to calculate residual fom glm on background or cueindex1/2
    cueindex=repmat([1 2 3 4],1,6); cueindex=repmat([1 2 1 2],1,6);
    
    for nn = size(site,1):-1:1,[aa,bb,cc]=glmfit([ones(1,size(cue_meanBK_tr(nn,:),2));squeeze(cue_meanBK_tr(nn,:))]',cue_mean_tr(nn,:)','normal');cueresid(nn,:)=cc.resid;end
    % for nn = 1:size(site,1),[aa,bb,cc]=glmfit([cueindex, cueindex2],squeeze(cue_mean_tr(nn,:)),'normal');cueresid(nn,:)=cc.resid;,end
    cueresid(isnan(cueresid))=0;
    %% decide which control we're using
    ccount=0;
    for control_flag = {'raw','bsub','bregress'}
        ccount = ccount+1;
        scatterindex = scatterindex+1;
        switch control_flag{1}
            case 'raw'
                cue_mean_tr2 = cue_mean_tr;
            case 'bsub'
                cue_mean_tr2 = cue_mean_tr-cue_meanBK_tr;
            case 'bregress'
                cue_mean_tr2 = cueresid;
        end
        figure;subplot(3,1,1),imagesc(corr(zscore(cue_mean_tr')')),subplot(3,1,2),imagesc(corr(zscore(cue_mean_tr'-cue_meanBK_tr')')),subplot(3,1,3),imagesc(corr(zscore(cueresid')'))

    iindex=[];for ii = [3,4,1,2],iindex = [iindex,ii,ii+4,ii+8,ii+12,ii+16,ii+20];end,
    
    figure(1111);xxx=corr(zscore(cue_mean_tr')');subplot(2,1,fi),imagesc(xxx(iindex,iindex))
        %% classifying on pseudoensembles
        
        
        clear x2 c2conf
        
        numsims = 1000;
        ncount=0;
        neuronvector=[5,10,25,50,75,100];
       % neuronvector=neuronvector(neuronvector<size(site,1));
        
        for numneurons=neuronvector
            ncount = ncount+1;
            ecount=0;
            for ee = numsims:-1:1
                
                ecount = ecount+1;
                ensembleorder=randperm(size(site,1),numneurons);
                %ensembleorder=randi(size(site,1),numneurons,1);
                
                
                testensemble = cue_mean_tr2(:,:);
                testensemble(isnan(testensemble))=0;
                
                bb = (testensemble(ensembleorder(1:numneurons),:)');
                numpcs=numneurons;
                %
                switch classifytype
                    case 'trial_blocks'
                        for trial = 6:-1:1
                            traintrials = setdiff(1:24,(1:4)+trial*4-4);
                            testtrials = (1:4)+trial*4-4;
                            
                            numpcs = numneurons;
                            
                            x2(ecount,:) = classify(bb(testtrials,1:numpcs),bb(traintrials,1:numpcs),cueindex(traintrials),'diagLinear');
                            
                            c2conf(ecount,1:4,1:4,trial)=zeros(4,4);
                            for ii = 1:size(x2,2)
                                c2conf(ecount,ii,x2(ecount,ii),trial)=...
                                    c2conf(ecount,ii,x2(ecount,ii),trial)+1;
                            end
                        end
                        
                    case 'trial_by_trial'
                        for trial = 24:-1:1
                            
                            traintrials = setdiff(1:24,trial);
                            testtrials = trial;
                            
                            numpcs = numneurons;
                            
                            x2(ecount,:) = classify(bb(testtrials,1:numpcs),bb(traintrials,1:numpcs),cueindex(traintrials),'diagLinear');
                            
                            c2conf(ecount,1:4,1:4,trial)=zeros(4,4);
                            
                            c2conf(ecount,1:4,1:4,trial)=zeros(4,4);
                            c2conf(ecount,cueindex(trial),x2(ecount,ii),trial)=...
                                c2conf(ecount,cueindex(trial),x2(ecount,ii),trial)+1;
                            
                        end
                end
            end
            classconf(ncount,ccount,fi,1:4,1:4,1:ecount)=shiftdim(squeeze(nansum(c2conf,4)),1);
        end
        
        %% classifying on the actual ensembles
        
        clear x2 c2conf
        
        ecount=0;
        for ee = unique(neuronEnsNum)'
            ecount = ecount+1;
            
            
            figure(100*ee);
            
            [aa, bb, cc]=princomp(zscore(cue_mean_tr2(neuronEnsNum==ee,:)'));
            hold on
            try
                subplot(2,3,scatterindex)
                title(strcat(['day ',num2str(fi),' animal ',num2str(ee),' condition ',control_flag{1}]))
                
                scatter(bb(:,1),bb(:,2),[],cueindex,'filled')
            catch
            end
            
            
            switch classifytype
                case 'trial_blocks'
                    for trial = 6:-1:1
                        
                        traintrials = setdiff(1:24,(1:4)+trial*4-4);
                        testtrials = (1:4)+trial*4-4;
                        
                        x2(ecount,:) = classify(cue_mean_tr2(neuronEnsNum==ee,testtrials)',cue_mean_tr2(neuronEnsNum==ee,traintrials)',cueindex(traintrials),'diagLinear');
                        
                        c2conf(ecount,1:4,1:4,trial)=zeros(4,4);
                        for ii = 1:size(x2,2)
                            c2conf(ecount,ii,x2(ecount,ii),trial)=...
                                c2conf(ecount,ii,x2(ecount,ii),trial)+1;
                        end
                    end
                    
                    if strcmp(control_flag{1},'raw')
                    try
                    [aa,bb,cc]=princomp(zscore(cue_mean_tr2(neuronEnsNum==ee,:)'));
                    
                    trial = 1;
                    
                    traintrials = setdiff(1:24,[]);
                    testtrials = (1:4)+trial*4-4;
                    
                    numpcs=2;
                    test = bb(testtrials,1:numpcs);
                    train = bb(traintrials,1:numpcs);
                    group = cueindex(traintrials);
                    
                    plot_lda_class(train,test,group,fi,ee)
                    catch
                    end
                    end
                    
                case 'trial_by_trial'
                    
                    for trial = 24:-1:1
                        
                        traintrials = setdiff(1:24,trial);
                        testtrials = trial;
                        
                        x2(ecount,:) = classify(cue_mean_tr2(neuronEnsNum==ee,testtrials)',cue_mean_tr2(neuronEnsNum==ee,traintrials)',cueindex(traintrials),'diagLinear');
                        
                        c2conf(ecount,1:4,1:4,trial)=zeros(4,4);
                        c2conf(ecount,cueindex(trial),x2(ecount,ii),trial)=...
                            c2conf(ecount,cueindex(trial),x2(ecount,ii),trial)+1;
                        
                    end
            end
        end
        classconfENS(fi,ccount,1:4,1:4,1:ecount)=shiftdim(squeeze(sum(c2conf,4)),1);
        ensnum(fi)={unique(neuronEnsNum)'};
        
        
        %% classifying on the actual neurons
        
        clear x2 c2conf
        
        ecount=0;
        for ee = 1:size(site,1)
            ecount = ecount+1;
            
            switch classifytype
                case 'trial_blocks'
                    
                    for trial = 6:-1:1
                        
                        traintrials = setdiff(1:24,(1:4)+trial*4-4);
                        testtrials = (1:4)+trial*4-4;
                        
                        x2(ecount,:) = classify(cue_mean_tr2(ee,testtrials)',cue_mean_tr2(ee,traintrials)',cueindex(traintrials),'diagLinear');
                        
                        c2conf(ecount,1:4,1:4,trial)=zeros(4,4);
                        for ii = 1:size(x2,2)
                            c2conf(ecount,ii,x2(ecount,ii),trial)=...
                                c2conf(ecount,ii,x2(ecount,ii),trial)+1;
                        end
                    end
                    
                case 'trial_by_trial'
                    
                    for trial = 24:-1:1
                        
                        traintrials = setdiff(1:24,trial);
                        testtrials = trial;
                        
                        x2(ecount,:) = classify(cue_mean_tr2(ee,testtrials)',cue_mean_tr2(ee,traintrials)',cueindex(traintrials),'diagLinear');
                        
                        c2conf(ecount,1:4,1:4,trial)=zeros(4,4);
                        c2conf(ecount,cueindex(trial),x2(ecount,ii),trial)=...
                            c2conf(ecount,cueindex(trial),x2(ecount,ii),trial)+1;
                        
                    end
                    
            end
        end
        
        classconfNRN(fi,ccount,1:4,1:4,1:ecount)=shiftdim(squeeze(sum(c2conf,4)),1);
        sizeSite(fi) = size(site,1);
        
    end
    
end

%% PLOTTING THE OUTPUT

clear betw* withinp* samep*
pltct =0;

for ctrl_num=3:-1:1
    for nn = 6:-1:1
        for ii = 2:-1:1,
            betweenpair(nn,ii,:)=nanmean([squeeze(classconf(nn,ctrl_num,ii,2,3,:)),squeeze(classconf(nn,ctrl_num,ii,3,2,:)),squeeze(classconf(nn,ctrl_num,ii,1,4,:)),squeeze(classconf(nn,ctrl_num,ii,4,1,:)),squeeze(classconf(nn,ctrl_num,ii,2,1,:)),squeeze(classconf(nn,ctrl_num,ii,1,2,:)),squeeze(classconf(nn,ctrl_num,ii,3,4,:)),squeeze(classconf(nn,ctrl_num,ii,4,3,:))],2);
        end,
    end
    for nn = 6:-1:1
        for ii = 2:-1:1,
            withinpair(nn,ii,:)=nanmean([squeeze(classconf(nn,ctrl_num,ii,1,3,:)),squeeze(classconf(nn,ctrl_num,ii,3,1,:)),squeeze(classconf(nn,ctrl_num,ii,2,4,:)),squeeze(classconf(nn,ctrl_num,ii,4,2,:))],2);
        end,
    end
    for nn = 6:-1:1
        for ii = 2:-1:1,
            samepair(nn,ii,:)=nanmean([squeeze(classconf(nn,ctrl_num,ii,1,1,:)),squeeze(classconf(nn,ctrl_num,ii,2,2,:)),squeeze(classconf(nn,ctrl_num,ii,3,3,:)),squeeze(classconf(nn,ctrl_num,ii,4,4,:))],2);
        end
    end
    clr = [.5 .5 .5; 0 0 0];
    figure(12145+1)
    for ii = 1:2
        subplot(3,3,1+(ctrl_num-1)*3)
        errorbar(neuronvector,nanmean(betweenpair(:,ii,:),3),nanstd(betweenpair(:,ii,:),[],3),'color',clr(ii,:)), hold on,ylim([0 6]),xlim([.5 100.5]),box off
        subplot(3,3,2+(ctrl_num-1)*3)
        errorbar(neuronvector,nanmean(withinpair(:,ii,:),3),nanstd(withinpair(:,ii,:),[],3),'color',clr(ii,:)), hold on,ylim([0 6]),xlim([.5 100.5]),box off,
        subplot(3,3,3+(ctrl_num-1)*3)
        errorbar(neuronvector,nanmean(samepair(:,ii,:),3),nanstd(samepair(:,ii,:),[],3),'color',clr(ii,:)), hold on,ylim([0 6]),xlim([.5 100.5]),box off,
    end
    
    histvec=-3.125:.25:3;
    pltct=pltct+2;
    
    for xx = 5:-1:1
        test_stat(xx,1,ctrl_num) = sum((squeeze(betweenpair(xx,2,:))'-squeeze(withinpair(xx,2,:))')>=0)/1000;
        test_stat(xx,2,ctrl_num) = sum((squeeze(withinpair(xx,1,:))'-squeeze(withinpair(xx,2,:))')>=0)/1000;
        
        figure(12341+xx);
        subplot(3,2,pltct-1),h=bar(histvec,histc(squeeze(withinpair(xx,1,:))'-squeeze(withinpair(xx,2,:))',histvec));shading flat,box off,xlim([-3.5 3.5]),ylim([0 270]),set(h,'barWidth',.8,'facecolor',[.2 .2 .2])
        hold on,text(.30,150,['p = ',num2str(test_stat(xx,2,ctrl_num))])
        hold on,text(0+.30,250,'day 1 > day 2')
        hold on,text(-3.50+.30,250,'day 2 > day 1')
        
        plot ([0 0],[0 250],'r--')
        
        subplot(3,2,pltct),h=bar(histvec,histc(squeeze(betweenpair(xx,2,:))'-squeeze(withinpair(xx,2,:))',histvec));shading flat,box off,xlim([-3.5 3.5]),ylim([0 270]),set(h,'barWidth',.8,'facecolor',[.2 .2 .2])
        hold on,text(.30,150,['p = ',num2str(test_stat(xx,1,ctrl_num))])
        hold on,text(0+.30,250,'unpaired cue >')
        hold on,text(-3.5+.30,250,'paired cue >')
        plot ([0 0],[0 250],'r--')
        
        
    end
end

for ctrl_num=3:-1:1
    clear betw* withinp* samep*
    
    for ii = 1:2,
        betweenpair(ii,:)=nansum([squeeze(classconfENS(ii,ctrl_num,2,3,:)),squeeze(classconfENS(ii,ctrl_num,3,2,:)),squeeze(classconfENS(ii,ctrl_num,1,4,:)),squeeze(classconfENS(ii,ctrl_num,4,1,:)),squeeze(classconfENS(ii,ctrl_num,2,1,:)),squeeze(classconfENS(ii,ctrl_num,1,2,:)),squeeze(classconfENS(ii,ctrl_num,3,4,:)),squeeze(classconfENS(ii,ctrl_num,4,3,:))],2)/8;
    end,
    
    for ii = 1:2,
        withinpair(ii,:)=nansum([squeeze(classconfENS(ii,ctrl_num,1,3,:)),squeeze(classconfENS(ii,ctrl_num,3,1,:)),squeeze(classconfENS(ii,ctrl_num,2,4,:)),squeeze(classconfENS(ii,ctrl_num,4,2,:))],2)/4;
    end
    for ii = 1:2,
        samepair(ii,:)=nansum([squeeze(classconfENS(ii,ctrl_num,1,1,:)),squeeze(classconfENS(ii,ctrl_num,2,2,:)),squeeze(classconfENS(ii,ctrl_num,3,3,:)),squeeze(classconfENS(ii,ctrl_num,4,4,:))],2)/4;
    end
    
    %
    % figure
    ens_bothdays = intersect(ensnum{1},ensnum{2});
    for ee= length(ens_bothdays):-1:1
        for ii = 2:-1:1
            within_ens(ii,ee) = withinpair(ii,ensnum{ii}==ens_bothdays(ee))/(4*6);
        end
    end
 
    
    betweenpair(betweenpair==0)=NaN;
    withinpair(withinpair==0)=NaN;
    samepair(samepair==0)=NaN;
    
    
    ens_bothdays = intersect(ensnum{1},ensnum{2});
    for ee= length(ens_bothdays):-1:1
        for ii = 2:-1:1
            between_ens(ii,ee) = betweenpair(ii,ensnum{ii}==ens_bothdays(ee))/(8*6);
        end
    end
    
    [hval(ctrl_num),pval(ctrl_num)]=ttest(betweenpair(2,:),withinpair(2,:));
    [hval2(ctrl_num),pval2(ctrl_num)]=ttest(withinpair(1,:),withinpair(2,:));
    
end

hval,pval,hval2,pval2


figure,
subplot(3,1,1)
errorbar(nanmean(betweenpair,2),nanstd(betweenpair,[],2)), hold on
bar(nanmean(betweenpair,2))

subplot(3,1,2)
errorbar(nanmean(withinpair,2),nanstd(withinpair,[],2)), hold on
bar(nanmean(withinpair,2))
subplot(3,1,3)
errorbar(nanmean(samepair,2),nanstd(samepair,[],2)), hold on,
bar(nanmean(samepair,2))


for ctrl_num = 3:-1:1
    clear betw* withinp* samep*
    for ii = 1:2,
        betweenpair(ii,:)=nansum([squeeze(classconfNRN(ii,ctrl_num,2,3,:)),squeeze(classconfNRN(ii,ctrl_num,3,2,:)),squeeze(classconfNRN(ii,ctrl_num,1,4,:)),squeeze(classconfNRN(ii,ctrl_num,4,1,:)),squeeze(classconfNRN(ii,ctrl_num,2,1,:)),squeeze(classconfNRN(ii,ctrl_num,1,2,:)),squeeze(classconfNRN(ii,ctrl_num,3,4,:)),squeeze(classconfNRN(ii,ctrl_num,4,3,:))],2)/8;
    end,
    
    for ii = 1:2,
        withinpair(ii,:)=nansum([squeeze(classconfNRN(ii,ctrl_num,1,3,:)),squeeze(classconfNRN(ii,ctrl_num,3,1,:)),squeeze(classconfNRN(ii,ctrl_num,2,4,:)),squeeze(classconfNRN(ii,ctrl_num,4,2,:))],2)/4;
    end
    for ii = 1:2,
        samepair(ii,:)=nansum([squeeze(classconfNRN(ii,ctrl_num,1,1,:)),squeeze(classconfNRN(ii,ctrl_num,2,2,:)),squeeze(classconfNRN(ii,ctrl_num,3,3,:)),squeeze(classconfNRN(ii,ctrl_num,4,4,:))],2)/4;
    end
    
    betweenpair(betweenpair==0)=NaN;
    withinpair(withinpair==0)=NaN;
    samepair(samepair==0)=NaN;
    
    figure,
    subplot(3,1,1)
    errorbar(nanmean(betweenpair,2),nanstd(betweenpair,[],2))
    subplot(3,1,2)
    errorbar(nanmean(withinpair,2),nanstd(withinpair,[],2))
    subplot(3,1,3)
    errorbar(nanmean(samepair,2),nanstd(samepair,[],2))
        [hval(ctrl_num),pval(ctrl_num)]=ttest(betweenpair(2,:),withinpair(2,:));
    [hval2(ctrl_num),pval2(ctrl_num)]=ttest(withinpair(1,:),withinpair(2,:));

end

hval,pval,hval2,pval2
