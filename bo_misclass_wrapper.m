function [classify_output, classify_output2] = bo_misclass_wrapper(neuron_index,...
    psth_AB_trials, psth_CD_trials, steps,...
    numperms, enssizes, max_A_trial, max_B_trial,window)

nindex = find(neuron_index);
ee=length(enssizes);
for enssize=enssizes(length(enssizes):-1:1)       
        for ii = numperms:-1:1
            %rp_index = randperm(length(nindex));
            rp_index = randi(length(nindex),1,enssize);

            thepsth2=psth_CD_trials(nindex(rp_index(1:enssize)),3:4,1:max_A_trial,:);
            thepsth1=psth_AB_trials(nindex(rp_index(1:enssize)),1:2,1:max_B_trial,:);
            [classify_output(ii,ee,:), classify_output2(ii,ee).blerg] = bo_psthmisclassify_constA(thepsth1,thepsth2,steps,window);
        end
        ee = ee-1;
end