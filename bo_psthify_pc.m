function [meanpsth, bsubpsth, AUC, steps, psth_AB_trials, psth_CD_trials,psth_AB_background, psth_CD_background] = bo_psthify_pc(site,neuronEnsNum,bin)


startime=bin; %analysis start (seconds)
endtime=80;   %analysis stop (seconds)
bstartime = bin;
bendtime =20;

steps=(startime:bin:endtime); %bins to use for PSTH & plots
cues = 1:4; %stimuli to include
ABtrials = 6;
CDtrials = 6; %# of trials to include

clear psth*
clear mean*
clear CUEZ
clear bsub*

AUC = zeros(size(site,1),size(site,2),length(startime:bin:endtime));

for i = [size(site,1):-1:1]
    l = 0;
    ens = neuronEnsNum(i);
    for j = cues
        
        if j<2.5
            for k = ABtrials:-1:1
                try
           
                    startime = bin+10;
                    endtime = 80;
                    
                    l=l+1;
                    psth_AB_trials(i,j,k,:) = [histc(site(i,j).trials(k).times,[startime:bin:endtime])];
                    psth_AB_background(i,j,k,:) = [histc(site(i,j).trials(k).times,[bstartime:bin:bendtime])];
                    psth_all_trials(i,l,:) = [histc(site(i,j).trials(k).times,[startime:bin:endtime])];
                    
                catch
                    psth_AB_trials(i,j,k,:) = NaN;
                    psth_AB_background(i,j,k,:) = NaN;
                    psth_all_trials(i,l,:) = NaN;
                    
                end
            end
        else
            
            j2=j-2;
            startime=bin; %analysis start (seconds)
            endtime=70;
            for k = CDtrials:-1:1
                try
                    l=l+1;
                    psth_all_trials(i,l,:) = [histc(site(i,j2).trials(k).times,[startime:bin:endtime])];
                    psth_CD_trials(i,j,k,:) = [histc(site(i,j2).trials(k).times,[startime:bin:endtime])];
                    psth_CD_background(i,j,k,:) = [histc(site(i,j2).trials(k).times,[bstartime:bin:bendtime])];
                    
                    
                    
                catch
                    
                    psth_CD_trials(i,j,k,:) = NaN;
                    psth_CD_background(i,j,k,:) = NaN;
                    psth_all_trials(i,l,:) = NaN;
                    
                    
                end
            end
        end
    end
    
    % calculate auROC for each stimulus
    for j = cues
        for t = size(psth_AB_trials,4):-1:1
            if j<2.5
                AUC(i,j,t) = auROC(reshape(squeeze(psth_AB_background(i,j,:,:)),1,[]),squeeze(psth_AB_trials(i,j,:,t)));
            else
                AUC(i,j,t) = auROC(reshape(squeeze(psth_CD_background(i,j,:,:)),1,[]),squeeze(psth_CD_trials(i,j,:,t)));
            end
            
        end
    end
end



%% calculate mean psth across stimuli
for i = size(site,1):-1:1
    l = 0;
    for j = cues
        if j<2.5
            meanpsth(i,j,:) =  squeeze(nanmean(psth_AB_trials(i,j,:,:),3));
        else
            meanpsth(i,j,:) =  squeeze(nanmean(psth_CD_trials(i,j,:,:),3));
            
        end
    end
end

%
AUC(isnan(AUC))=.5;

backgroundHz = squeeze(mean(meanpsth(:,:,find(steps==bstartime,1,'first'):find(steps==bendtime,1,'first')),3));

for nn = 1:size(site,1)
    for cc = 1:size(site,2)
        
        bsubpsth(nn,cc,:)=meanpsth(nn,cc,:)-backgroundHz(nn,cc);
    end
end

steps = steps-40;

