
function [probeResults, trialprobe] = bo_behavior_probe(pokein,pokeout,pcflag)

for hh = 1:size(pokein,1)
    for ii = 1:4 %FOR ALL 4 CUES
        
        jjct = 0;
        if ii<2.5
            
            cueon =40;
            cueoff=50;
            baseon =35;
            baseoff=40;
            
            if pcflag
            cueon =50;
            cueoff=60;
            end
            
            for jj = 1:min([length(pokein(hh,ii).trials),length(pokeout(hh,ii).trials),6])
                
                jjct=jjct+1;
                
                %%% POKE DURATION PER TRIAL
                
                
                pokeduration{hh,ii}(jj) = 0;
                pokebaseline{hh,ii}(jj) = 0;
                pokecounter{hh,ii}(jj) = 0;
                
                for kk = 1:length(find(pokein(hh,ii).trials(jj).times<80))
                    
                    outindex = find(pokeout(hh,ii).trials(jj).times>pokein(hh,ii).trials(jj).times(kk),1,'first');
                    
                    if ~isempty(outindex)
                      
                        outtime = pokeout(hh,ii).trials(jj).times(outindex);
                        
                        if kk<length(find(pokein(hh,ii).trials(jj).times<70))
                            if outtime>pokein(hh,ii).trials(jj).times(kk+1)%if there are 2 pokeins in a row
                                outtime = pokein(hh,ii).trials(jj).times(kk)+.0001;
                            end
                        end
                        
                        inptime = pokein(hh,ii).trials(jj).times(kk);
                        %% FOR CUE DURATION CALCULATION
                        if inptime<cueon && outtime>cueoff
                            pokeduration{hh,ii}(jj) = cueoff-cueon;
                        elseif inptime<cueon && outtime>cueon && outtime<cueoff
                            pokeduration{hh,ii}(jj) = pokeduration{hh,ii}(jj)+ outtime-cueon;
                        elseif inptime>cueon && inptime<cueoff && outtime>cueoff
                            pokecounter{hh,ii}(jj) = pokecounter{hh,ii}(jj)+1;
                            pokeduration{hh,ii}(jj) = pokeduration{hh,ii}(jj)+ cueoff - inptime;
                        elseif inptime>cueon && inptime<cueoff && outtime<cueoff && outtime>cueon
                            pokeduration{hh,ii}(jj) = pokeduration{hh,ii}(jj)+outtime-inptime;
                            pokecounter{hh,ii}(jj) = pokecounter{hh,ii}(jj)+1;
                        end
                        
                        %% FOR CUE BASELINE CALCULATION
                        if inptime<baseon && outtime>baseoff
                            pokebaseline{hh,ii}(jj) = cueoff-cueon;
                        elseif inptime<baseon && outtime>baseon && outtime<baseoff
                            pokebaseline{hh,ii}(jj) = pokebaseline{hh,ii}(jj)+ outtime-baseon;
                        elseif inptime>baseon && inptime<baseoff && outtime>baseoff
                            pokebaseline{hh,ii}(jj) = pokebaseline{hh,ii}(jj)+ baseoff - inptime;
                        elseif inptime>baseon && inptime<baseoff && outtime<baseoff && outtime>baseon
                            pokebaseline{hh,ii}(jj) = pokebaseline{hh,ii}(jj)+outtime-inptime;
                        end
                        
                    end
                end
            end
            
        elseif ii>2.5
            ii2 = ii;
            
            if pcflag
            cueon =40;
            cueoff=50;
            ii2 = ii-2;
            end
            
            for jj = 1:min([length(pokein(hh,ii2).trials),length(pokeout(hh,ii2).trials),6])
                
                jjct=jjct+1;
                
                %%% POKE DURATION PER TRIAL
                
                
                pokeduration{hh,ii}(jj) = 0;
                pokebaseline{hh,ii}(jj) = 0;
                pokecounter{hh,ii}(jj) = 0;
                
                for kk = 1:length(find(pokein(hh,ii2).trials(jj).times<80))
                    
                    outindex = find(pokeout(hh,ii2).trials(jj).times>pokein(hh,ii2).trials(jj).times(kk),1,'first');
                    
                    if ~isempty(outindex)
                        
                        outtime = pokeout(hh,ii2).trials(jj).times(outindex);
                        
                        if kk<length(find(pokein(hh,ii2).trials(jj).times<70))
                            if outtime>pokein(hh,ii2).trials(jj).times(kk+1)%if there are 2 pokeins in a row
                                outtime = pokein(hh,ii2).trials(jj).times(kk)+.0001;
                            end
                        end
                        
                        inptime = pokein(hh,ii2).trials(jj).times(kk);
                        %% FOR CUE DURATION CALCULATION
                        if inptime<cueon && outtime>cueoff
                            pokeduration{hh,ii}(jj) = cueoff-cueon;
                        elseif inptime<cueon && outtime>cueon && outtime<cueoff
                            pokeduration{hh,ii}(jj) = pokeduration{hh,ii}(jj)+ outtime-cueon;
                        elseif inptime>cueon && inptime<cueoff && outtime>cueoff
                            pokecounter{hh,ii}(jj) = pokecounter{hh,ii}(jj)+1;
                            pokeduration{hh,ii}(jj) = pokeduration{hh,ii}(jj)+ cueoff - inptime;
                        elseif inptime>cueon && inptime<cueoff && outtime<cueoff && outtime>cueon
                            pokeduration{hh,ii}(jj) = pokeduration{hh,ii}(jj)+outtime-inptime;
                            pokecounter{hh,ii}(jj) = pokecounter{hh,ii}(jj)+1;
                        end
                        
                        %% FOR CUE BASELINE CALCULATION
                        if inptime<baseon && outtime>baseoff
                            pokebaseline{hh,ii}(jj) = cueoff-cueon;
                        elseif inptime<baseon && outtime>baseon && outtime<baseoff
                            pokebaseline{hh,ii}(jj) = pokebaseline{hh,ii}(jj)+ outtime-baseon;
                        elseif inptime>baseon && inptime<baseoff && outtime>baseoff
                            pokebaseline{hh,ii}(jj) = pokebaseline{hh,ii}(jj)+ baseoff - inptime;
                        elseif inptime>baseon && inptime<baseoff && outtime<baseoff && outtime>baseon
                            pokebaseline{hh,ii}(jj) = pokebaseline{hh,ii}(jj)+outtime-inptime;
                        end
                        
                    end
                end
            end
        end
    end
    
end


for jj = 1:6
for hh = 1:size(pokein,1)
    try
        B(hh,1) = nanmean(pokebaseline{hh,1}(1:3));
        D(hh,1) = nanmean(pokebaseline{hh,2}(1:3));
    end
    try 
        A(hh,1,jj) = (pokebaseline{hh,3}(jj));
        C(hh,1,jj) = (pokebaseline{hh,4}(jj));
        
    end
    try
        B(hh,2) = nanmean(pokeduration{hh,1}(1:3));
        D(hh,2) = nanmean(pokeduration{hh,2}(1:3));
    end
    try 
        A(hh,2,jj) = (pokeduration{hh,3}(jj));
        C(hh,2,jj) = (pokeduration{hh,4}(jj));
    end
end
end

for hh = 1:size(pokein,1)
    for ii = 1:size(pokein,2)
        pcount(hh,ii) = sum(pokecounter{hh,ii}(1:end)>0);
    end
end
zz=cueoff-cueon;
for tr=1:6

    probeBaseline = [B(:,1), D(:,1), A(:,1,tr), C(:,1,tr)];
    probeDuration = [B(:,2), D(:,2), A(:,2,tr), C(:,2,tr)];
    probeResults = probeDuration-probeBaseline;

%  

figure(111)
subplot(1,3,1),bar(nanmean(probeResults(:,1:4)/zz)),hold on
errorbar(1:4,squeeze(nanmean(probeResults(:,1:4)/zz,1)),squeeze(nanstd(probeResults(:,1:4)/zz,1))/sqrt(22)*1)
%[xx yy]=ttest(probeDuration(:,3)-probeDuration(:,4));
trialprobe(:,3:4,tr) = probeResults(:,3:4);
ylabel('fraction of cue in food cup')
ylim([0 0.7])
end


probeDuration = [B(:,2), D(:,2), nanmean(A(:,2,:),3), nanmean(C(:,2,:),3)];
probeBasline = [B(:,1), D(:,1), nanmean(A(:,1,:),3), nanmean(C(:,1,:),3)];
probeResults = probeDuration;

[a b c]=anova1(probeResults)

%% within subjects error bars
for tt = 1:size(trialprobe,3)
         grand_average = nanmean(nanmean(trialprobe(:,[3 4],tt),1),2);

for ii = 1:size(trialprobe,1)

     individual_average = nanmean(trialprobe(ii,[3 4],tt),2);
     datatoplot(ii,[1 2],tt) = trialprobe(ii,[3 4],tt)-individual_average+grand_average;
end
end
 
 
 figure(111)
 subplot(1,3,2)
 errorbar(1:6,squeeze(nanmean(datatoplot(:,1,:),1))/zz,squeeze(nanstd(datatoplot(:,1,:),1)/zz)/sqrt(22)*1),hold on
 errorbar(1:6,squeeze(nanmean(datatoplot(:,2,:),1))/zz,squeeze(nanstd(datatoplot(:,2,:),1)/zz)/sqrt(22)*1)

 
xlim([.5 6.5])
xlabel('trial #')
ylabel('fraction of cue in food cup')
[signtest_p signtest_h]=signtest(nanmean(trialprobe(:,3,:),3),nanmean(trialprobe(:,4,:),3))
[ttest_h ttest_p x ttest_statAC]=ttest(nanmean(trialprobe(:,3,:),3),nanmean(trialprobe(:,4,:),3))
[ttest_h ttest_p x ttest_statBD]=ttest(nanmean(trialprobe(:,1,:),3),nanmean(trialprobe(:,2,:),3))
box off
ylim([0 0.7])


subplot(1,3,3)
bar([-.5:.025:.499],histc(nanmean(trialprobe(:,3,:),3)/zz-nanmean(trialprobe(:,4,:),3)/zz,[-.5:.025:-0.0001,0.0001:.025:.5]))
xlabel('fraction responding, A minus C')
ylabel('number of animals')
xlim([-.5 .5])

anoGroup = [];
group1var= [];group2var= [];

for ii = 1:6
   anoGroup = [anoGroup trialprobe(:,3:4,ii)'];
   group1var = [group1var repmat(ii,2,length(trialprobe(:,3,ii))')];
   group2var = [group2var repmat(3:4,length(trialprobe(:,3,ii))',1)'];
end



[a b c] = anovan(anoGroup(:),{group2var(:)},'model','full');

