function [meanpsth, bsubpsth, AUC, steps, psth_AB_trials, psth_CD_trials,psth_AB_background, psth_CD_background] = bo_psthify_probe(site,neuronEnsNum,bin)


startime=bin; %analysis start (seconds)
endtime=80;   %analysis stop (seconds)
bstartime = 10;
bendtime =40;

steps=(startime:bin:endtime); %bins to use for PSTH & plots
win = 10; % analysis window for each method below (# of BINS)

cues = 1:4; %stimuli to include
ABtrials = 3;
CDtrials = 2; %# of trials to include

clear psth*
clear mean*
clear CUEZ
clear bsub*

AUC = zeros(size(site,1),size(site,2),length(startime:bin:endtime));

%% bin spikes into PSTH
for i = [size(site,1):-1:1]
    l = 0;
    ens = neuronEnsNum(i);
    for j = cues
        
        if j<2.5
            for k = ABtrials:-1:1
                
%          
         
                try
                l=l+1;
                psth_AB_trials(i,j,k,:) = [histc(site(i,j).trials(k).times,[startime:bin:endtime])];
                psth_AB_background(i,j,k,:) = [histc(site(i,j).trials(k).times,[bstartime:bin:bendtime])];
                psth_all_trials(i,l,:) = [histc(site(i,j).trials(k).times,[startime:bin:endtime])];
                catch
                end
            end
            
        else
            for k = CDtrials:-1:1
                l=l+1;
                try
                psth_all_trials(i,l,:) = [histc(site(i,j).trials(k).times,[startime:bin:endtime])];
                psth_CD_trials(i,j,k,:) = [histc(site(i,j).trials(k).times,[startime:bin:endtime])];
                psth_CD_background(i,j,k,:) = [histc(site(i,j).trials(k).times,[bstartime:bin:bendtime])];
                %
                catch
                end

            end
        end
    end
    
    % calculate auROC for each stimulus
    for j = cues
        for t = size(psth_AB_trials,4):-1:1
            if j<2.5
                                 AUC(i,j,t) = auROC(reshape(squeeze(psth_AB_background(i,j,1:3,:)),1,[]),squeeze(psth_AB_trials(i,j,1:3,t)));
            else
                                 AUC(i,j,t) = auROC(reshape(squeeze(psth_CD_background(i,j,1:end,:)),1,[]),squeeze(psth_CD_trials(i,j,1:end,t)));
            end
            
        end
    end
end

AUC(isnan(AUC))=0.5;

%% calculate mean psth across stimuli
for i = size(site,1):-1:1
    l = 0;
    for j = cues
        if j<2.5
            meanpsth(i,j,:) =  squeeze(nanmean(psth_AB_trials(i,j,:,:),3));
            %meanPokepsth(i,j,:,:) =  squeeze(nanmean(psth_AB_poke(i,j,:,:),3));
            %             meanDippsth(i,j,:,:) =  squeeze(nanmean(psth_B_dipper(i,1,:,:,:),3));
        else
            meanpsth(i,j,:) =  squeeze(nanmean(psth_CD_trials(i,j,:,:),3));
            %             meanPokepsth(i,j,:,:) =  squeeze(nanmean(psth_CD_poke(i,j,:,:),3));
            
        end
    end
end

backgroundHz = squeeze(mean(meanpsth(:,:,find(steps==bstartime,1,'first'):find(steps==bendtime,1,'first')),3));

for nn = 1:size(site,1)
    for cc = 1:size(site,2)
        
        bsubpsth(nn,cc,:)=meanpsth(nn,cc,:)-mean(backgroundHz(nn,cc),2);
        
    end
end

steps = steps-40;

