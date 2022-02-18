% ccepDTI02_optimizeDetector
% determination of optimal setting for n1-detector for sEEG data

% author: Dorien van Blooijs & Susanne Jelsma
% date: October 2021

% load ECoG data, split into stimulation trials, re-reference the data,
% load scoring, calculate kappa, 
% SUMMARY
%% set paths
% set umcuEpi_CCEP_DTI/matlab in your directory and run this section

clc
clear
cfg.folderinput = 'chronic_ECoG'; % from which folder would you like to load ECoGs?
myDataPath = setLocalDataPath(cfg);

%% patient characteristics
%instead of these two lines search automaticaly for the visual scored patients
    %sub_label = input('Patient number (RESPXXXX) (select multiple patients by separating each with a comma): ','s'); %(select multiple patients by separating each with a comma)
    %cfg.sub_label = strsplit(sub_label,{', ',','});

if isfield(myDataPath,'CCEPpath2')
    files2 = dir(myDataPath.CCEPpath2);
else
    error('myDataPath.CCEPpath2 does not exist. Make sure you add this folder for a second observer in personalDataPath.m');
end

files1 = dir(myDataPath.CCEPpath);
files1(1:2) = []; files2(1:2) = [];
sub_label1 = cell(1,size(files1,1));
for subj=1:size(files1,1)
    sub_label1(1,subj)= cellstr(erase(files1(subj).name,'sub-'));
end

sub_label2 = cell(1,size(files2,1));
for subj=1:size(files2,1)
    sub_label2(1,subj)= cellstr(erase(files2(subj).name,'sub-'));
end
sub_label = intersect(sub_label1,sub_label2);

i = 1;
for subj= 1:size(sub_label,2)
    x = input(sprintf('Load subject %s? (y/n): ',char(sub_label(subj))),'s');       
    if strcmp(x,'y') 
        cfg.sub_label(1,i) = sub_label(subj) ;
        i=i+1;
    end
end

clear sub_label sub_label1 sub_label2 files1 files2 i x subj

cfg = selectPatients(cfg, myDataPath);


%% load ECoGs 

dataBase = load_ECoGdata(myDataPath,cfg);

disp('All ECoGs are loaded')

%% load visual scoring results

dataBase = load_visualScores(dataBase,myDataPath,cfg);

%% preprocessing CCEP in ECoG

% clear cfg

% preprocessing step
cfg.dir = 'no'; % if you want to take negative/positive stimulation into account
cfg.amp = 'no'; % if you want to take stimulation current into account

% select epochs and average
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

dataBase = preprocess_ECoG(dataBase, cfg);

disp('All ECoGs are preprocessed')

%% rereference data
% find 10 signals in the same trial with lowest variance and not being a
% bad channel and use that as reference signal to increase the
% signal-to-noise ratio

cfg.reref = input('Do you want to rereference the data with an average of the 10 signals with the lowest variance? [y/n]: ','s'); % ('yes' = re-reference, 'no' = no re-reference)

for subj = 1:size(dataBase,2)
    for run = 1:size(dataBase(subj).metadata,2)
        
        if strcmp(cfg.reref,'y')
            
            dataBase = rerefCCEPdata(dataBase,subj,run,cfg);

        else
            
            dataBase(subj).metadata(run).cc_epoch_sorted_reref = dataBase(subj).metadata(run).cc_epoch_sorted;
        end
        
        % take mean of these (not) re-referenced signals
        dataBase(subj).metadata(run).cc_epoch_sorted_reref_avg = squeeze(mean(dataBase(subj).metadata(run).cc_epoch_sorted_reref,2,'omitnan')); 
    
        fprintf('...%s %s has been re-referenced... \n',...
            dataBase(subj).sub_label,dataBase(subj).metadata(run).run_label)

    end
end

if strcmp(cfg.reref,'y')
    disp('Data is rereferenced and all trials of one stimulus pair are averaged.')
else
    disp('Data is not re-referenced and all trials of one stimulus pair are averaged.')
end


%% runs samenvoegen
for subj= 1:size(dataBase,2)
    if size(cfg.run_label{subj},2) > 1
       dataBase(subj).metadata_runs = dataBase(subj).metadata; % kopieer alles naar metadata_runs
       data = dataBase(subj).metadata(1); % laat alle info van run 1 de nieuwe metadata zijn
       dataBase(subj).metadata = data;
            for run = 2:size(cfg.run_label{subj},2)
                scored = dataBase(subj).metadata_runs(run).ccep_VS1.checked;
                scored2 = dataBase(subj).metadata_runs(run).ccep_VS2.checked; 
                data_sorted = dataBase(subj).metadata_runs(run).cc_epoch_sorted;
                data_sorted_reref = dataBase(subj).metadata_runs(run).cc_epoch_sorted_reref;
                data_sorted_avg = dataBase(subj).metadata_runs(run).cc_epoch_sorted_avg; 
                data_sorted_reref_avg = dataBase(subj).metadata_runs(run).cc_epoch_sorted_reref_avg;
                stimsets = dataBase(subj).metadata_runs(run).cc_stimsets;
                stimchannels = dataBase(subj).metadata_runs(run).cc_stimchans;
                    
    
                dataBase(subj).metadata.ccep_VS1.checked = [dataBase(subj).metadata.ccep_VS1.checked scored];
                dataBase(subj).metadata.ccep_VS2.checked = [dataBase(subj).metadata.ccep_VS2.checked scored2];
                dataBase(subj).metadata.cc_epoch_sorted_avg = [dataBase(subj).metadata.cc_epoch_sorted_avg data_sorted_avg];
                dataBase(subj).metadata.cc_epoch_sorted_reref_avg = [dataBase(subj).metadata.cc_epoch_sorted_reref_avg data_sorted_reref_avg];
                dataBase(subj).metadata.cc_epoch_sorted = cat(3,dataBase(subj).metadata.cc_epoch_sorted, data_sorted);
                dataBase(subj).metadata.cc_epoch_sorted_reref = cat(3,dataBase(subj).metadata.cc_epoch_sorted_reref, data_sorted_reref);
                dataBase(subj).metadata.cc_stimsets = [dataBase(subj).metadata.cc_stimsets;stimsets];
                dataBase(subj).metadata.cc_stimchans = [dataBase(subj).metadata.cc_stimchans;stimchannels];
                %voeg hier nog alle andere data die je gebruikt die per run
                %verschilt toe!!! vooral voor de rescore
            end
    end
end

disp('runs merged')
clear scored scored2 data data_sorted_reref_avg data_sorted_avg data_sorted data_sorted_reref run subj stimsets stimchannels 
%% kappa score

% pre-allocation
subjs = size(dataBase,2);
TP = NaN(subjs,1); % true positives
FN = NaN(subjs,1); % false negatives
FP = NaN(subjs,1); % false positives
TN = NaN(subjs,1); % true negatives   
kappa_score = NaN(subjs,1); % Cohen's kappa

for subj = 1:size(dataBase,2)
    scored = dataBase(subj).metadata.ccep_VS1.checked;
    scored2 = dataBase(subj).metadata.ccep_VS2.checked;
    TP(subj) = numel(find( scored2 == 1 & scored == 1));
    FN(subj) = numel(find( scored2 == 1 & scored == 0));
    FP(subj) = numel(find( scored2 == 0 & scored == 1));
    TN(subj) = numel(find( scored2 == 0 & scored == 0));
    kappa_score(subj) = (2*(TP(subj)*TN(subj)-FN(subj)*FP(subj)))...
        /((TP(subj)+FP(subj))*(FP(subj)+TN(subj))+(TP(subj)+FN(subj))*(FN(subj)+TN(subj))) ;
end
% with the external function kappa
kappa_score2 = NaN(subjs,1); % Cohen's kappa via function

for subj = 1:size(dataBase,2)
    scored = dataBase(subj).metadata.ccep_VS1.checked;
    scored2 = dataBase(subj).metadata.ccep_VS2.checked;
    TP(subj) = numel(find( scored2 == 1 & scored == 1));
    FN(subj) = numel(find( scored2 == 1 & scored == 0));
    FP(subj) = numel(find( scored2 == 0 & scored == 1));
    TN(subj) = numel(find( scored2 == 0 & scored == 0));
    confusion_matrix = [TP(subj),FP(subj);FN(subj),TN(subj)];
    kappa_score2(subj) = kappa(confusion_matrix);
    clear('confusion_matrix');
end
clear subjs subj scored2 scored TP FN FP TN
%% merge scores of two scorers
% only the scored CCEPs we both agreed on are officially visual scored
for subj = 1:size(dataBase,2)
    scored = dataBase(subj).metadata.ccep_VS1.checked;
    scored2 = dataBase(subj).metadata.ccep_VS2.checked;       
    scored_both = scored + scored2;
    scored_both(scored_both == 1) = 0; % when there is disagreement, value is 1, this should be regarded as no CCEP
    scored_both(scored_both == 2) = 1; % when there is agreement, value is 2, this should be regarded as a CCEP
    dataBase(subj).metadata.visual_scored = scored_both;
end

clear scored_both scored scored2 subj
disp('scores merged')
%% rescore CCEP with no consensus
% in een functie/ander script 
% % subj=3;
% % run=1;
% % 
% % filefolder = fullfile(myDataPath.CCEPpath,dataBase(subj).sub_label,dataBase(subj).ses_label,dataBase(subj).metadata(run).run_label);
% % if ~exist(filefolder,'dir')
% %     mkdir(filefolder)
% % end
% % 
% % scored = dataBase(subj).metadata.ccep_VS1.checked;
% % scored2 = dataBase(subj).metadata.ccep_VS2.checked;
% % scored_both = scored + scored2;
% % [RE_chan,RE_stimp] = find( scored_both == 1);
% % 
% % cfg.n1Detected = 'n'; % --> if n1 peaks are detected
% % rescored = dataBase(subj).metadata.visual_scored;
% % 
% % 
% % for i = 1:length(RE_chan)
% %     stimp = RE_stimp(i);
% %     chan = RE_chan(i);
% %     fs = dataBase(subj).metadata(run).ccep_header.Fs;
% %     tt = -cfg.epoch_prestim+1/fs:1/fs:cfg.epoch_length-cfg.epoch_prestim;
% %     
% %     H = plot_ccep(dataBase,subj,run,cfg,tt,chan,stimp);
% %     
% %     perc = i / length(RE_chan) *100;
% %     x = input(sprintf('%2.1f %% --- stimpair = %s-%s chan = %s --- Is this a CCEP? (y/n): ',...
% %     perc,dataBase(subj).metadata(run).cc_stimchans{stimp,:},dataBase(subj).metadata(run).ch{chan}),'s');
% %     % rescore the CCEPs
% %     
% %     if strcmp(x,'y') 
% %     rescored(chan,stimp) = 1 ;
% %     else
% %     rescored(chan,stimp) = 0 ;
% %     end
% %     close(H)
% %     clear('x')
% %     % save also till which stimpair visual N1s are checked.
% %     cfg.checkUntilStimp = stimp;
% %     
% %     filename = [dataBase(subj).sub_label,'_',dataBase(subj).ses_label,'_',dataBase(subj).metadata(run).task_label,'_',dataBase(subj).metadata(run).run_label,'_N1sREChecked.mat'];
% %     
% %     % save file during scoring in case of error
% %     save(fullfile(filefolder,filename),'rescored');    
% % end
% % dataBase(subj).metadata.visual_scored = rescored;   
% % fprintf('Rescore of %s completed\n',dataBase(subj).sub_label)

% %%
% clear filename filefolder rescored scored scored2 scored_both stimp subj tt run RE_chan RE_stimp x i H fs chan perc
%% load rescored data for use again
dataPath = myDataPath.CCEPpath; 
for subj = 1:size(dataBase,2)
    sub_label = ['sub-' cfg.sub_label{subj}];
    ses_label = cfg.ses_label{subj};
    task_label = cfg.task_label{subj};
    run_label = cfg.run_label{subj}{1};
    DRE = dir(fullfile(dataPath,sub_label,ses_label,run_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label,'_N1sREChecked.mat'])); % rescored data
    if size(DRE,1) == 0
            error('%s does not exist',dir(fullfile(dataPath,sub_label,ses_label,run_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label,'_N1sREChecked.mat'])));
    end 
    dataNameRE = fullfile(DRE(1).folder, DRE(1).name);
    dataRE = load(dataNameRE);
    dataBase(subj).metadata.visual_scored = dataRE.rescored;
end
disp('rescored data loaded')
clear dataNameRE dataRE DRE run_label ses_label sub_label task_label subj dataPath
%% optimize detector
% in een ander script want run je niet elke keer?
% for cECoG data, these are the best parameters:
% % cfg.amplitude_thresh = 2.6;
% % cfg.n1_peak_range = 100;
% % cfg.minSD = 50;
% % cfg.sel = 20;
% % 
% % % range of parameters tested:
% % amplitude_tresh_range = 0.5:0.1:1 ;
% % minSD_range = 90:1:110 ;
% % sel_range = 0:1:20 ;
% % 
% % n=1;
% % % pre-allocation
% % subjs = size(dataBase,2);
% % combs = size(minSD_range,2)*size(sel_range,2)*size(amplitude_tresh_range,2);
% % TP = NaN(combs,subjs); % true positives
% % FN = NaN(combs,subjs); % false negatives
% % FP = NaN(combs,subjs); % false positives
% % TN = NaN(combs,subjs); % true negatives
% % combination = NaN(combs,3); % combination of parameters. dit geeft een error
% % 
% % for amplTh = amplitude_tresh_range 
% % for sd = minSD_range 
% % for sl = sel_range
% %     cfg.amplitude_thresh = amplTh;
% %     cfg.minSD = sd;
% %     cfg.sel = sl;
% %     dataBase = detect_n1peak_ccep(dataBase, cfg); 
% %     combination(n,:) = [amplTh sd sl]; % combination(n= number combination, parameter) parameter 1= amplitude_thresh; 2=minSD; 3=sel
% %     for subj = 1:size(dataBase,2)
% %         detected = dataBase(subj).metadata.ccep.n1_peak_sample; 
% %         detected(~isnan(detected)) = 1;
% %         detected(isnan(detected)) = 0;
% %         scored = dataBase(subj).metadata.visual_scored;
% % 
% %         TP(n,subj) = numel(find( scored == 1 & detected == 1)); 
% %         FN(n,subj) = numel(find( scored == 1 & detected == 0));
% %         FP(n,subj) = numel(find( scored == 0 & detected == 1));
% %         TN(n,subj) = numel(find( scored == 0 & detected == 0));
% %     end
% %     n=n+1;
% %     fprintf('...Combination: Amplitude Threshold=%g, min SD of baseline=%g, Sel=%g, has been run...\n',...
% %         amplTh,sd,sl)
% % end
% % end
% % end
%clear n subjs combs sel_range minSD_range amplitude_tresh_range detected scored sd sl amplTh


%% 2.4 without sel 0-10 
% combination=combination(in_combs,:);
% FN=FN(in_combs,:);
% FP=FP(in_combs,:);
% TN=TN(in_combs,:);
% TP=TP(in_combs,:);
%% calculate and visualize performance for each setting
% analyze the performances per subject
subjs = size(dataBase,2);
combs = size(combination,1);
nr_visual_scored = NaN(subjs,1); % check if nr FP/FN/TN/TP is same as nr visual scored
spec_subj = NaN(combs,subjs); % specificity
sens_subj = NaN(combs,subjs); % sensitivity
ppv_subj = NaN(combs,subjs); % positive predictive value
npv_subj = NaN(combs,subjs); % negative predicitive value
d_roc_subj = NaN(combs,subjs); % distance to top left corner of ROC 
d_prc_subj = NaN(combs,subjs); % distance to top left corner of PRC (precision-recall curve
F_score_subj = NaN(combs, subjs); % F-score (wss hetzelfde als d_roc)
F2_score_subj = NaN(combs, subjs); % F2-score  sensitivity to be twice as important as positive predictive value 
beta=2;

for subj = 1:size(dataBase,2)
   nr_visual_scored(subj) = size(dataBase(subj).metadata.visual_scored,1)*size(dataBase(subj).metadata.visual_scored,2);
   for n = 1:length(combination)
        spec_subj(n,subj) = TN(n,subj)/(TN(n,subj)+FP(n,subj));
        sens_subj(n,subj) = TP(n,subj)/(TP(n,subj)+FN(n,subj));
        ppv_subj(n,subj) = TP(n,subj)/(TP(n,subj)+FP(n,subj));
        npv_subj(n,subj) = TN(n,subj)/(TN(n,subj)+FN(n,subj));
        d_roc_subj(n,subj) = sqrt((1-sens_subj(n,subj))^2+ (1-spec_subj(n,subj))^2);
        d_prc_subj(n,subj) = sqrt((1-sens_subj(n,subj))^2+ (1-ppv_subj(n,subj))^2);
        F_score_subj(n,subj) = TP(n,subj)/(TP(n,subj)+0.5*(FP(n,subj)+FN(n,subj)));
        F2_score_subj(n,subj) = ((1+beta^2)*TP(n,subj))/((1+beta^2)*TP(n,subj)+(beta^2)*FN(n,subj)+FP(n,subj));
   end
end

% combine the performances of each subj and run.
TP_all = NaN(combs,1);
FN_all = NaN(combs,1);
FP_all = NaN(combs,1);
TN_all = NaN(combs,1);
spec = NaN(combs,1); % specificity
sens = NaN(combs,1); % sensitivity
ppv = NaN(combs,1); % positive predictive value
npv = NaN(combs,1); % negative predicitive value
d_roc = NaN(combs,1); % distance to top left corner of ROC 
d_prc = NaN(combs,1); % distance to top left corner of PRC (precision-recall curve)
F_score = NaN(combs,1); % F-score (wss hetzelfde als d_roc)
F2_score = NaN(combs, 1); % F2-score  sensitivity to be twice as important as positive predictive value 

for n = 1:length(combination)
    TP_all(n,1) = sum(TP(n,:), 'omitnan');
    FN_all(n,1) = sum(FN(n,:),'omitnan');
    FP_all(n,1) = sum(FP(n,:),'omitnan');
    TN_all(n,1) = sum(TN(n,:),'omitnan');
    if TP_all(n) + FN_all(n) + FP_all(n) + TN_all(n) == sum(nr_visual_scored)
        spec(n,1) = TN_all(n)/(TN_all(n)+FP_all(n));
        sens(n,1) = TP_all(n)/(TP_all(n)+FN_all(n));
        ppv(n,1) = TP_all(n)/(TP_all(n)+FP_all(n));
        npv(n,1) = TN_all(n)/(TN_all(n)+FN_all(n));
        d_roc(n,1) = sqrt((1-sens(n))^2+ (1-spec(n))^2);
        d_prc(n,1) = sqrt((1-sens(n))^2+ (1-ppv(n))^2);  
        F_score(n,1) = TP_all(n)/(TP_all(n)+0.5*(FP_all(n)+FN_all(n)));
        F2_score(n,1) = ((1+beta^2)*TP(n))/((1+beta^2)*TP(n)+(beta^2)*FN(n)+FP(n));
    else
        error('Number of visual scored CCEPs is not the same as sum of all detected');
    end
end
clear subjs subj n beta
%% best combination for this range of parameters
[value, best_comb] = min(d_prc);
best_amplTh = combination(best_comb,1);
best_minSD = combination(best_comb,2);
best_sel = combination(best_comb,3);
%%
% plot ROC and PRC
f1 = figure(1);
scatter(1-spec,sens, 30, [0 0.4470 0.7410], '.')
hold on
scatter(1-spec(best_comb), sens(best_comb), 50, [0.8500 0.3250 0.0980],'filled')
title('ROC curve')
xlabel(' 1-specificity')
ylabel(' sensitivity')
xlim([0,1])
ylim([0,1])
ticks = 0:0.1:1;
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f1,'-dpng', 'ROC_d_prc_2_5','-r300')
saveas(f1,'ROC_d_prc_2_5.m')

f2 = figure(2);
scatter(sens,ppv, 30, [0 0.4470 0.7410], '.')
hold on
scatter(sens(best_comb), ppv(best_comb), 50, [0.8500 0.3250 0.0980], 'filled')
title('Precision-recall curve')
xlabel('sensitivity')
ylabel('positive predictive value')
xlim([0,1])
ylim([0,1])
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f2,'-dpng', 'PRC__d_prc_2_5','-r300')
saveas(f2,'PRC_d_prc_2_5.m')
%% more sensitive algorithm
% best combination for this range of parameters
[value2, best_comb2] = min(d_roc);
best_amplTh2 = combination(best_comb2,1);
best_minSD2 = combination(best_comb2,2);
best_sel2 = combination(best_comb2,3);
%%
% plot ROC and PRC
f1 = figure(1);
scatter(1-spec,sens, 30, [0 0.4470 0.7410], '.')
hold on
scatter(1-spec(best_comb2), sens(best_comb2), 50, [0.8500 0.3250 0.0980],'filled')
title('ROC curve')
xlabel(' 1-specificity')
ylabel(' sensitivity')
xlim([0,1])
ylim([0,1])
ticks = 0:0.1:1;
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f1,'-dpng', 'ROC_d_roc_2_5','-r300')
saveas(f1,'ROC_d_roc_2_5.m')

f2 = figure(2);
scatter(sens,ppv, 30, [0 0.4470 0.7410], '.')
hold on
scatter(sens(best_comb2), ppv(best_comb2), 50, [0.8500 0.3250 0.0980], 'filled')
title('Precision-recall curve')
xlabel('sensitivity')
ylabel('positive predictive value')
xlim([0,1])
ylim([0,1])
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f2,'-dpng', 'PRC__d_roc_2_5','-r300')
saveas(f2,'PRC_d_roc_2_5.m')
%% algorithm based on f2score
% best combination for this range of parameters
[value3, best_comb3] = max(F2_score);
best_amplTh3 = combination(best_comb3,1);
best_minSD3 = combination(best_comb3,2);
best_sel3 = combination(best_comb3,3);
%%
% plot ROC and PRC
f1 = figure(1);
scatter(1-spec,sens, 30, [0 0.4470 0.7410], '.')
hold on
scatter(1-spec(best_comb3), sens(best_comb3), 50, [0.8500 0.3250 0.0980],'filled')
title('ROC curve')
xlabel(' 1-specificity')
ylabel(' sensitivity')
xlim([0,1])
ylim([0,1])
ticks = 0:0.1:1;
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f1,'-dpng', 'ROC_F2_score_2_5','-r300')
saveas(f1,'ROC_F2_score_2_5.m')

f2 = figure(2);
scatter(sens,ppv, 30, [0 0.4470 0.7410], '.')
hold on
scatter(sens(best_comb3), ppv(best_comb3), 50, [0.8500 0.3250 0.0980], 'filled')
title('Precision-recall curve')
xlabel('sensitivity')
ylabel('positive predictive value')
xlim([0,1])
ylim([0,1])
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f2,'-dpng', 'PRC_F2_score_2_5','-r300')
saveas(f2,'PRC_F2_score_2_5.m')

%% algorithm based on fscore
% best combination for this range of parameters
[value4, best_comb4] = max(F_score);
best_amplTh4 = combination(best_comb4,1);
best_minSD4 = combination(best_comb4,2);
best_sel4 = combination(best_comb4,3);
%%
% plot ROC and PRC
f1 = figure(1);
scatter(1-spec,sens, 30, [0 0.4470 0.7410], '.')
hold on
scatter(1-spec(best_comb4), sens(best_comb4), 50, [0.8500 0.3250 0.0980],'filled')
title('ROC curve')
xlabel(' 1-specificity')
ylabel(' sensitivity')
xlim([0,1])
ylim([0,1])
ticks = 0:0.1:1;
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f1,'-dpng', 'ROC_F1_score_2_5','-r300')
saveas(f1,'ROC_F1_score_2_5.m')

f2 = figure(2);
scatter(sens,ppv, 30, [0 0.4470 0.7410], '.')
hold on
scatter(sens(best_comb4), ppv(best_comb4), 50, [0.8500 0.3250 0.0980], 'filled')
title('Precision-recall curve')
xlabel('sensitivity')
ylabel('positive predictive value')
xlim([0,1])
ylim([0,1])
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f2,'-dpng', 'PRC_F1_score_2_5','-r300')
saveas(f2,'PRC_F1_score_2_5.m')

%% all performance parameters in one plot all subj
% plot ROC and PRC
f1 = figure(1);
scatter(1-spec,sens, 30, [0 0.4470 0.7410], '.')
hold on
scatter(1-spec(best_comb2), sens(best_comb2), 50, [0.8500 0.3250 0.0980],'filled')
scatter(1-spec(best_comb), sens(best_comb), 100, [0.9290 0.6940 0.1250],'filled')
scatter(1-spec(best_comb4), sens(best_comb4), 50, [0.4940 0.1840 0.5560],'filled')
scatter(1-spec(best_comb3), sens(best_comb3), 50, [0.4660 0.6740 0.1880],'filled')
title('ROC curve')
xlabel(' 1-specificity')
ylabel(' sensitivity')
legend('range of parameters','d-roc','d-prc','F1-score','F2-score','Location','southeast')
xlim([0,1])
ylim([0,1])
ticks = 0:0.1:1;
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f1,'-dpng', 'ROC_subjall_2_1','-r300')
saveas(f1,'ROC_subjall_2_1.m')

f2 = figure(2);
scatter(sens,ppv, 30, [0 0.4470 0.7410], '.')
hold on
scatter(sens(best_comb2), ppv(best_comb2), 50, [0.8500 0.3250 0.0980], 'filled')
scatter(sens(best_comb), ppv(best_comb), 100, [0.9290 0.6940 0.1250], 'filled')
scatter(sens(best_comb4), ppv(best_comb4), 50, [0.4940 0.1840 0.5560], 'filled')
scatter(sens(best_comb3), ppv(best_comb3), 50, [0.4660 0.6740 0.1880], 'filled')
title('Precision-recall curve')
xlabel('sensitivity')
ylabel('positive predictive value')
legend('range of parameters','d-roc','d-prc','F1-score','F2-score','Location','southwest')
xlim([0,1])
ylim([0,1])
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f2,'-dpng', 'PRC_subjall_2_1','-r300')
saveas(f2,'PRC_subjall_2_1.m')

%% all performance parameters in one plot per subj
% plot ROC and PRC
subj=3;
f1 = figure(1);
scatter(1-spec_subj(:,subj),sens_subj(:,subj), 30, [0 0.4470 0.7410], '.')
hold on
scatter(1-spec_subj(best_comb2,subj), sens_subj(best_comb2,subj), 100, [0.8500 0.3250 0.0980],'filled')
scatter(1-spec_subj(best_comb,subj), sens_subj(best_comb,subj), 100, [0.9290 0.6940 0.1250],'filled')
scatter(1-spec_subj(best_comb4,subj), sens_subj(best_comb4,subj), 50, [0.4940 0.1840 0.5560],'filled')
scatter(1-spec_subj(best_comb3,subj), sens_subj(best_comb3,subj), 50, [0.4660 0.6740 0.1880],'filled')
title('ROC curve')
xlabel(' 1-specificity')
ylabel(' sensitivity')
legend('range of parameters subject 3','d-roc','d-prc','F1-score','F2-score','Location','southeast')
xlim([0,1])
ylim([0,1])
ticks = 0:0.1:1;
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f1,'-dpng', 'ROC_subj3_2_3','-r300')
saveas(f1,'ROC_subj3_2_3.m')

f2 = figure(2);
scatter(sens_subj(:,subj),ppv_subj(:,subj), 30, [0 0.4470 0.7410], '.')
hold on
scatter(sens_subj(best_comb2,subj), ppv_subj(best_comb2,subj), 100, [0.8500 0.3250 0.0980], 'filled')
scatter(sens_subj(best_comb,subj), ppv_subj(best_comb,subj), 100, [0.9290 0.6940 0.1250], 'filled')
scatter(sens_subj(best_comb4,subj), ppv_subj(best_comb4,subj), 50, [0.4940 0.1840 0.5560], 'filled')
scatter(sens_subj(best_comb3,subj), ppv_subj(best_comb3,subj), 50, [0.4660 0.6740 0.1880], 'filled')
title('Precision-recall curve')
xlabel('sensitivity')
ylabel('positive predictive value')
legend('range of parameters subject 3','d-roc','d-prc','F1-score','F2-score','Location','southwest')
xlim([0,1])
ylim([0,1])
set(gca, 'YTick',ticks, 'XTick', ticks);
box on
print(f2,'-dpng', 'PRC_subj3_2_3','-r300')
saveas(f2,'PRC_subj3_2_3.m')

%% overview of performance values of best combination for this range of parameters
best = best_comb2;
pf_value = { 'Specificity';'Sensitivity';'Positive predicitve value';'Negative predicitve value';'Distance to corner ROC (receiver operating characteristic) curve';'Distance to corner precision-recall curve (PRC)';'F1 score';'F2 score'};
pf_overal = round([spec(best,:); sens(best,:); ppv(best,:);npv(best,:);d_roc(best,:);d_prc(best,:);F_score(best,:);F2_score(best,:)].*100);
pf_subj = round([spec_subj(best,:); sens_subj(best,:);ppv_subj(best,:);npv_subj(best,:);d_roc_subj(best,:); d_prc_subj(best,:);F_score_subj(best,:);F2_score_subj(best,:)].*100);
pf_table = table(pf_value,pf_overal,pf_subj(:,1),pf_subj(:,2),pf_subj(:,3),'VariableNames',{'Performance','Overall','Subject 1','Subject 2','Subject 3'});



%% range of parameters 2.4 without sel 0-10

%[amplTh sd sl]
combs_ones_sel=find(combination(:,3)>=10);
combs_zeros_sel=find(combination(:,3)<10);
in_combs = combs_ones_sel;
out_combs = combs_zeros_sel;

%%
% plot ROC and PRC
f1 = figure(1);
scatter(1-spec(out_combs),sens(out_combs), 30, [0.3010 0.7450 0.9330], 'o')
hold on
scatter(1-spec(in_combs),sens(in_combs), 30, [0 0.4470 0.7410], '.')
scatter(1-spec(best_comb2), sens(best_comb2), 50, [0.8500 0.3250 0.0980],'filled')
%scatter(1-spec(best_comb), sens(best_comb), 100, [0.9290 0.6940 0.1250],'filled')
%scatter(1-spec(best_comb4), sens(best_comb4), 50, [0.4940 0.1840 0.5560],'filled')
%scatter(1-spec(best_comb3), sens(best_comb3), 50, [0.4660 0.6740 0.1880],'filled')
scatter(1-spec(2175), sens(2175), 50, [0.4660 0.6740 0.1880],'filled')
title('ROC curve','Fontsize',14)
xlabel(' 1-specificity','Fontsize',13)
ylabel(' sensitivity','Fontsize',13)
legend('range of parameters sel<10','range of parameters sel>l0','d-roc sel<10','d-roc sel>10','Location','southeast','Fontsize',12) %'d-prc','F1-score','F2-score'
%xlim([0,1])
%ylim([0,1])
ticks = 0:0.05:1;
set(gca, 'YTick',ticks, 'XTick', ticks,'Fontsize',12);
set(gcf, 'Position', get(0, 'Screensize'));
box on
print(f1,'-dpng', 'ROC_with_without_sel','-r300')
saveas(f1,'ROC_with_without_sel.m')

f2 = figure(2);
scatter(sens(out_combs),ppv(out_combs), 30, [0.3010 0.7450 0.9330], 'o')
hold on
scatter(sens(in_combs),ppv(in_combs), 30, [0 0.4470 0.7410], '.') % opactity dat je het door elkaar heen ziet
scatter(sens(best_comb2), ppv(best_comb2), 50, [0.8500 0.3250 0.0980], 'filled')
% scatter(sens(best_comb), ppv(best_comb), 100, [0.9290 0.6940 0.1250], 'filled')
% scatter(sens(best_comb4), ppv(best_comb4), 50, [0.4940 0.1840 0.5560], 'filled')
% scatter(sens(best_comb3), ppv(best_comb3), 50, [0.4660 0.6740 0.1880], 'filled')
scatter(sens(2175), ppv(2175), 50, [0.4660 0.6740 0.1880], 'filled')
title('Precision-recall curve','Fontsize',14)
xlabel('sensitivity','Fontsize',13)
ylabel('positive predictive value','Fontsize',13)
legend('range of parameters sel<10','range of parameters sel>l0','d-roc sel<10','d-roc sel>10','Location','southwest','Fontsize',12) 
%xlim([0,1])
%ylim([0,1])
set(gca, 'YTick',ticks, 'XTick', ticks,'Fontsize',12);
set(gcf, 'Position', get(0, 'Screensize'));
box on
print(f2,'-dpng', 'PRC_with_without_sel','-r300')
saveas(f2,'PRC_with_without_sel.m')


%% determine optimal settings


%% analyseer uitkomsten
 %droc evt nog F2 score?
% zet dit in command window: detect_n1peak_ccep
cfg.amplitude_thresh = combination(3508,1); % best_comb2 verander dit ook als je de F2 score veranderd!
cfg.minSD = combination(3508,2); % 2175 3508
cfg.sel = combination(3508,3);
cfg.n1_peak_range = 100;
cfg.n1Detected = 'n'; % --> if n1 peaks are detected
%% analyse peaks
subj=1;%voeg subject nog to aan detect_n1peak_ccep_analyse
chan= 13 ;
stimp= 1;
dataBase = detect_n1peak_ccep_analyse(dataBase, cfg,chan,stimp);
%%
dataBase = detect_n1peak_ccep(dataBase, cfg);

%%
subj=3;
run=1;
cfg.n1Detected = 'y';
fs = dataBase(subj).metadata(run).ccep_header.Fs;
tt = -cfg.epoch_prestim+1/fs:1/fs:cfg.epoch_length-cfg.epoch_prestim;

scored = dataBase(subj).metadata.visual_scored;
detected = dataBase(subj).metadata.ccep.n1_peak_sample;
detected(~isnan(detected)) = 1;
detected(isnan(detected)) = 0;
[FP_chan,FP_stimp]= find( scored == 0 & detected == 1);
%% analyse FP
for i = 100:length(FP_chan)
    stimp = FP_stimp(i);
    chan = FP_chan(i);
    H = plot_ccep(dataBase,subj,run,cfg,tt,chan,stimp);
    perc = i / length(FP_chan) *100;
    x = input(sprintf('%2.1f %% --- stimpair = %s-%s chan = %s --- save this false positive? (y/n): ',...
    perc,dataBase(subj).metadata(run).cc_stimchans{stimp,:},dataBase(subj).metadata(run).ch{chan}),'s');
end

%%
cfg.n1Detected = 'n';
[FN_chan,FN_stimp]= find( scored == 1 & detected == 0);
%% analyse FN
for i = 120:length(FN_chan)
    stimp = FN_stimp(i);
    chan = FN_chan(i);
    H = plot_ccep(dataBase,subj,run,cfg,tt,chan,stimp);
    perc = i / length(FN_chan) *100;
    x = input(sprintf('%2.1f %% --- stimpair = %s-%s chan = %s --- save this false negative? (y/n): ',...
    perc,dataBase(subj).metadata(run).cc_stimchans{stimp,:},dataBase(subj).metadata(run).ch{chan}),'s');
end
