function [output, output2] = psthmisclassify(thepsth1,thepsth2,steps)%%% develop a neural classifier based on a 'leave one out' approach on a
%%% per neuron basis with the input of psth_AB_trials


%% first, reshape data into matrices


win =6;

for nn = size(thepsth1,1):-1:1
    cue_by_trial_counter=0;
    for cue = size(thepsth1,2):-1:1
        for tr = size(thepsth1,3):-1:1
            
            cue_by_trial_counter = cue_by_trial_counter+1;
            
            all_psth1(cue_by_trial_counter,nn,:)=thepsth1(nn,cue,tr,:);
                        
            cuegroup_test(cue_by_trial_counter) = cue;
            trialgroup_test(cue_by_trial_counter) = tr;
            
        end
    end
end


for nn = size(thepsth1,1):-1:1
    cue_by_trial_counter=0;
    for cue = size(thepsth2,2):-1:1
        for tr = size(thepsth2,3):-1:1
            
            cue_by_trial_counter = cue_by_trial_counter+1;
            
            all_psth2(cue_by_trial_counter,nn,:)=thepsth2(nn,cue,tr,:);
            
            cuegroup_train(cue_by_trial_counter) = cue;
            trialgroup_train(cue_by_trial_counter) = tr;
            
        end
    end
end

clear numcorrect
for tbin = ceil(win/2)+1:size(thepsth1,4)-ceil(win/2)-1
    dim = 1;
    if max(dim>size(thepsth1,1))
        dim = 1:size(thepsth1,1);
    end
    
    [aa bb]=princomp(([...
        nanmean(all_psth1(:,:,tbin-ceil(win/2):tbin+ceil(win/2)),3);...
        nanmean(all_psth2(:,:,logical((steps>0).*(steps<10))),3)]')');
    
    
    test_pca = bb(1:length(cuegroup_test),dim);

    train_pca = bb(length(cuegroup_test)+1:length(cuegroup_test)+length(cuegroup_train),dim);

    %for the current neuron
    test_trial=1;
    
    try [xx]=classify(test_pca(:,:),train_pca(:,:),cuegroup_train(:),'diagLinear');
        %numcorrect(tbin,test_trial) = sum(xx==cuegroup(:))/size(thepsth1,2);
        %numcorrect(tbin,test_trial) = sum(xx==cuegroup([4 5 6 10 11 12])')/size(thepsth1,2);
        %numcorrect(tbin,test_trial) = sum(xx==cuegroup([1 2 3 7 8 9])')/size(thepsth1,2);
        numcorrect(tbin,test_trial) = sum(xx==cuegroup_test(:))/size(thepsth1,2);
        numcorrect1(tbin,test_trial) = sum(xx(cuegroup_test==1)==1);
        numcorrect2(tbin,test_trial) = sum(xx(cuegroup_test==2)==2);

    catch ME
        numcorrect(tbin,test_trial) = NaN;
    end
end



%figure(555)
%hold on

%plot(steps(1:end-4),medfilt1((numcorrect'),1)/size(thepsth1,3))

output = (numcorrect');
output2 = [numcorrect1';numcorrect2'];
