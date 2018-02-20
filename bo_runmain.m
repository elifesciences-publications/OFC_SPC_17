 %% OFC SPC MASTER SCRIPT


clear all

%% hard coded values -- sorting threshold // bin size // etc.
sortthreshold = 3; % remove neurons w/ deviation less than 3x noise
bin = .25;

%% choose the data you're going to load, and load it

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for preconditioning - un-comment these two files 
% and run script across both days worth of data, once

%filename = ('PC1_wo25.mat');pcflag = 1;
%filename2 = ('PC2_wo25.mat');pcflag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for probe data - un commment this file 

filename = ('probe_resort_wodrift_woobviousMUA3.mat');pcflag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for conditioning data, un-comment these files one at a time
% and run script for each day

%filename = ('C1_cueAligned_minus25-malfunction.mat'); pcflag=2;
%filename = ('C2_cueAligned_minus25_all21.mat'); pcflag=2;
%filename = ('C3_cueAligned_minus25_allbut17.mat'); pcflag=2;
%filename = ('C4_cueAligned_minus25_all21.mat'); pcflag=2;
%filename = ('C5_cueAligned_minus25_all21.mat'); pcflag=2;
%filename = ('C6_cueAligned_minus25_missing9.mat'); pcflag=2;

[site1, neuronEnsNum, pokein, pokeout, beh_ens] = load_spc(filename);

if exist('filename2','var')
    [site2, neuronEnsNum2, pokein, pokeout, beh_ens] = load_spc(filename2);
    neuronEnsNum2=neuronEnsNum2; 
    site = [site1;site2]; neuronEnsNum =[neuronEnsNum;neuronEnsNum2];
    neuronDayNum = zeros(size(neuronEnsNum));
    neuronDayNum(1:length(neuronEnsNum))=1;
    neuronDayNum(length(neuronEnsNum)-length(neuronEnsNum2):length(neuronEnsNum))=2;
else
    site=site1;
    neuronDayNum = zeros(size(neuronEnsNum));
    neuronDayNum(1:length(neuronEnsNum))=1;
end

%%  first remove neurons that don't exceed a signal-to-noise threshold
for ii = size(site,1):-1:1, sortstat(ii)=site(ii,1).stats.sig2noise;end
site = site(sortstat>sortthreshold,:);
neuronEnsNum = neuronEnsNum(sortstat>sortthreshold,:);
neuronDayNum = neuronDayNum(sortstat>sortthreshold,:);


 site=site(neuronEnsNum~=25,:);
 neuronEnsNum=neuronEnsNum(neuronEnsNum~=25,:);
 beh_ens=beh_ens(beh_ens~=25);

 pokein=pokein(beh_ens~=25,:);
  pokeout=pokeout(beh_ens~=25,:);


%% analyze behavior

% first extract, plot, and test relevent behavioral variables
[probeDuration, trialprobe] = bo_behavior_probe(pokein, pokeout, pcflag);

% with this data, if it's a probe day, assign grouping for good or poor
% performing animals

if pcflag % if it's a preconditioning session, do nothing.
else % if it's a probeday
    
    % pull out behavioral response for each animal
    % broad agreement between behavioral response measures, but multiple
    % measures were explored to ensure robustness 
    
    reg_y1 = nanmean(trialprobe(:,3,:),3)-nanmean(trialprobe(:,4,:),3); %% mean time A>C
    
    % find the neurons from good or bad behaved animals
    clear nindex*
    for ii = length(neuronEnsNum):-1:1,nindex1(ii)=sum(neuronEnsNum(ii)==beh_ens(reg_y1>0));end
    for ii = length(neuronEnsNum):-1:1,nindex2(ii)=sum(neuronEnsNum(ii)==beh_ens(reg_y1<=0));end

    nindex1 = logical(nindex1);  % the animals with 'good' behavior
    nindex2 = logical(nindex2);    % the animals with 'poor' behavior
    
    [beh_num,beh_order]=sort(reg_y1); % from low to high responding to A
end

%% psth_ify neural data and calculate basic stats

% make simple psths from neural data - for 3 measures - raw
% trialxtrial, raw mean psth, basline normalized psth, AUC normalized
% psth for all cues -- and provide a time index -- steps -- that gives
% cue-onset at time 0 (usually +/-40s from cue onset)


if pcflag == 1 % preconditioning
    
    nindex = 1:size(site,1);
    [meanpsth, bsubpsth, AUC, steps, psth_AB_trials,psth_CD_trials,psth_AB_background,psth_CD_background]= bo_psthify_pc(site,neuronEnsNum,bin);
    [numnrn_all] = bo_ttest_cueresponse(psth_AB_background,psth_AB_trials,psth_CD_background,psth_CD_trials,nindex,steps);
    [numnrn_all] = bo_ttest_cueresponse_inhib(psth_AB_background,psth_AB_trials,psth_CD_background,psth_CD_trials,nindex,steps);

elseif pcflag == 2 %conditioning
    
    [meanpsth, bsubpsth, AUC, steps, psth_AB_trials,psth_AB_background]= bo_psthify_cond(site,neuronEnsNum,bin);
    nindex = 1:size(site,1);
    [numnrn_all] = bo_ttest_cueresponse2(psth_AB_background,psth_AB_trials,nindex,steps);
    
else % probe test
    
    [meanpsth, bsubpsth, AUC, steps, psth_AB_trials,psth_CD_trials,psth_AB_background,psth_CD_background]= bo_psthify_probe(site,neuronEnsNum,bin);
    
    % find and plot fraction of responsive neurons for a session
    
    nindex = 1:size(site,1);
    
    [numnrn_good,nrnsgood] = bo_ttest_cueresponse3(psth_AB_background,psth_AB_trials,psth_CD_background,psth_CD_trials,nindex1,steps);
    [numnrn_bad,nrnsbad] = bo_ttest_cueresponse3(psth_AB_background,psth_AB_trials,psth_CD_background,psth_CD_trials,nindex2,steps);
    [numnrn_all,nrnsall] = bo_ttest_cueresponse3(psth_AB_background,psth_AB_trials,psth_CD_background,psth_CD_trials,nindex,steps);
    
   
end
%% further analysis - misclassification and correlation for probe and preconditioning data
% for either preconditioning or probe data, call relevent scripts to plot
% panels for figures 2 and 3 or figures 5 and 6, respectively



if pcflag == 1  %% if data is from a preconditioning day
    
    nindex = ones(size(site,2),1); neural_resp = AUC; %#ok<*NASGU,*PREALL>
    bo_psth_line_plot
    
    nindex =  1:size(site,1);  neural_resp = meanpsth;
    bo_rho_resample
    
    bo_corr_earlylate % for multibin correlation suppliment to figure 3
    
    bo_classify_PC
    
elseif pcflag == 2
else
    
    
    nindex = ones(size(site,2),1); neural_resp = AUC;
    bo_psth_line_plot
    
    nindex = nindex1; neural_resp = meanpsth;
    bo_rho_resample
    
    nindex = nindex2; neural_resp = meanpsth;
    bo_rho_resample
    
    bo_probe_classification
    
end