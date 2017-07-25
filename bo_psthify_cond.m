function [meanpsth, bsubpsth, AUC, steps, psth_AB_trials, psth_AB_background] = bo_psthify_probe(site,neuronEnsNum,bin)


startime=bin; %analysis start (seconds)
endtime=80;   %analysis stop (seconds)
bstartime = bin;
bendtime =30;

steps=(startime:bin:endtime); %bins to use for PSTH & plots
win = 10; % analysis window for each method below (# of BINS)

cues = 1:2; %stimuli to include
ABtrials = 6;
CDtrials = 6; %# of trials to include

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
                
%                 if j == 1
%                     ztimes = [1 4.5 8];
%                     for zz = 1:3
%                         psth_B_dipper(i,j,k,zz,:) = [histc(site(i,j).trials(k).times,[startime:bin:endtime]-ztimes(zz))];
%                     end
%                 end
                
                % ptimes = pokeIns(ens,j).trials(k).times;
                
                % if ~isempty(ptimes)
                % psth_AB_firstcuepoke(i,j,k,:) =  [histc(site(i,j).trials(k).times,[ptimes(1)-5:bin:ptimes(1)+5])];
                
                %    for pp = 1:length(ptimes)
                %       ptime = ptimes(pp);
                %       psth_AB_pokes(pp,:) =  [histc(site(i,j).trials(k).times,[ptime-5:bin:ptime+5])];
                %  end
                % psth_AB_poke(i,j,k,:)=squeeze(mean(psth_AB_pokes,1));
                % end
                try
                l=l+1;
                psth_AB_trials(i,j,k,:) = [histc(site(i,j).trials(k).times,[startime:bin:endtime])];
                psth_AB_background(i,j,k,:) = [histc(site(i,j).trials(k).times,[bstartime:bin:bendtime])];
                psth_all_trials(i,l,:) = [histc(site(i,j).trials(k).times,[startime:bin:endtime])];
                catch
                end
            end
            
        else
  
        end
    end
    
    % calculate auROC for each stimulus
    for j = cues
        for t = size(psth_AB_trials,4):-1:1
            if j<2.5
                                 AUC(i,j,t) = auROC(reshape(squeeze(psth_AB_background(i,j,:,:)),1,[]),squeeze(psth_AB_trials(i,j,:,t)));
            else
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
            
        end
    end
end

backgroundHz = squeeze(mean(meanpsth(:,:,find(steps==bstartime,1,'first'):find(steps==bendtime,1,'first')),3));

for nn = 1:size(site,1)
    for cc = 1:size(site,2)
        
        bsubpsth(nn,cc,:)=meanpsth(nn,cc,:)-backgroundHz(nn,cc);
    end
end

steps = steps-40;

