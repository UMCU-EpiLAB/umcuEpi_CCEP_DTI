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
%instead of these two lines search automaticly for the visual scored patients
    %sub_label = input('Patient number (RESPXXXX) (select multiple patients by separating each with a comma): ','s'); %(select multiple patients by separating each with a comma)
    %cfg.sub_label = strsplit(sub_label,{', ',','});

files1 = dir(myDataPath.CCEPpath); files2 = dir(myDataPath.CCEPpath2);
files1(1:2) = []; files2(1:2) = [];
sub_label1 = cell(1,size(files1,1));
for subj=1:size(files1,1)
sub_label1(1,subj)= cellstr(erase(files1(subj).name,'sub-'));
end
sub_label2 = cell(1,size(files1,1));
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

cfg = selectPatients(cfg, myDataPath);

%% load visual scoring results

dataBase_visualScores = load_visualScores(myDataPath,cfg);

%% load ECoGs 

dataBase = load_ECoGdata(myDataPath,cfg);

disp('All ECoGs are loaded')

%% preprocessing CCEP in ECoG

clear cfg

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

%% kappa score

% pre-allocation
TP = NaN(1); % true positives
FN = NaN(1); % false negatives
FP = NaN(1); % false positives
TN = NaN(1); % true negatives   
kappa_score = NaN(1); % Cohen's kappa

    for subj = 1:size(dataBase_visualScores,2)
        for run = 1:size(dataBase_visualScores(subj).metadata,2)
            scored = dataBase_visualScores(subj).metadata(run).ccep.checked;
            scored2 = dataBase_visualScores2(subj).metadata(run).ccep.checked;
            TP(subj,run) = numel(find( scored2 == 1 & scored == 1));
            FN(subj,run) = numel(find( scored2 == 1 & scored == 0));
            FP(subj,run) = numel(find( scored2 == 0 & scored == 1));
            TN(subj,run) = numel(find( scored2 == 0 & scored == 0));
            kappa_score(subj,run) = (2*(TP(subj,run)*TN(subj,run)-FN(subj,run)*FP(subj,run)))...
                /((TP(subj,run)+FP(subj,run))*(FP(subj,run)+TN(subj,run))+(TP(subj,run)+FN(subj,run))*(FN(subj,run)+TN(subj,run))) ;
        end
    end

% with the external function kappa
kappa_score2 = NaN(1); % Cohen's kappa via function

    for subj = 1:size(dataBase_visualScores,2)
        for run = 1:size(dataBase_visualScores(subj).metadata,2)
            scored = dataBase_visualScores(subj).metadata(run).ccep.checked;
            scored2 = dataBase_visualScores2(subj).metadata(run).ccep.checked;
            TP(subj,run) = numel(find( scored2 == 1 & scored == 1));
            FN(subj,run) = numel(find( scored2 == 1 & scored == 0));
            FP(subj,run) = numel(find( scored2 == 0 & scored == 1));
            TN(subj,run) = numel(find( scored2 == 0 & scored == 0));
            confusion_matrix = [TP(subj,run),FP(subj,run);FN(subj,run),TN(subj,run)];
            kappa_score2(subj,run) = kappa(confusion_matrix);
            clear('confusion_matrix');
        end
    end

%% merge scores of two scorers
% only the scored CCEPs we both agreed on are officially visual scored
    for subj = 1:size(dataBase_visualScores,2)
        for run = 1:size(dataBase_visualScores(subj).metadata,2)
            scored = dataBase_visualScores(subj).metadata(run).ccep.checked;
            scored2 = dataBase_visualScores(subj).metadata(run).ccep2.checked;       
            scored_both = scored + scored2;
            scored_both(scored_both == 1) = 0; % when there is disagreement, value is 1, this should be regarded as no CCEP
            scored_both(scored_both == 2) = 1; % when there is agreement, value is 2, this should be regarded as a CCEP
            dataBase_visualScores(subj).metadata(run).visual_scored = scored_both;
        end
    end

clear scored_both scored scored2 dataBase_visualScores2
%% optimize detector
% for cECoG data, these are the best parameters:
cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;
cfg.minSD = 50;
cfg.sel = 20;

n=1;
% pre-allocation
TP = NaN(1); % true positives
FN = NaN(1); % false negatives
FP = NaN(1); % false positives
TN = NaN(1); % true negatives
%combination = NaN(1); % combination of parameters. dit geeft een error

for sd = 11:1:13 
for sl = 10:1:12 
for amplTh = 3.2:0.1:3.3 
    cfg.amplitude_thresh = amplTh;
    cfg.minSD = sd;
    cfg.sel = sl;
    dataBase = detect_n1peak_ccep(dataBase, cfg);
    combination(n,:) = [amplTh sd sl]; % combination(n= number combination, parameter) parameter 1= amplitude_thresh; 2=minSD; 3=sel
    for subj = 1:size(dataBase,2)
        for run = 1:size(dataBase(subj).metadata,2)
            detected = dataBase(subj).metadata(run).ccep.n1_peak_sample; 
            detected(~isnan(detected)) = 1;
            detected(isnan(detected)) = 0;
            scored = dataBase_visualScores(subj).metadata(run).visual_scored;

            TP(subj,run,n) = numel(find( scored == 1 & detected == 1)); 
            FN(subj,run,n) = numel(find( scored == 1 & detected == 0));
            FP(subj,run,n) = numel(find( scored == 0 & detected == 1));
            TN(subj,run,n) = numel(find( scored == 0 & detected == 0));
            
        end
    end
    n=n+1;
    fprintf('...Combination: Amplitude Threshold=%s, min SD of baseline=%s, Sel=%s, has been run...\n',...
        amplTh,sd,sl)
end
end
end

%% calculate and visualize performance for each setting

% combine the performances of each subj and run.
TP_all = NaN(1);
FN_all = NaN(1);
FP_all = NaN(1);
TN_all = NaN(1);
spec = NaN(1); % specificity
sens = NaN(1); % sensitivity
ppv = NaN(1); % positive predictive value
npv = NaN(1); % negative predicitive value
d_roc = NaN(1); % distance to top left corner of ROC 
d_prc = NaN(1); % distance to top left corner of PRC (precision-recall curve)
F_score = NaN(1); % F-score (wss hetzelfde als d_roc)

        for n = 1:length(combination)
        TP_all(n,1) = sum(nonzeros(TP(:,:,n)));
        FN_all(n,1) = sum(nonzeros(FN(:,:,n)));
        FP_all(n,1) = sum(nonzeros(FP(:,:,n)));
        TN_all(n,1) = sum(nonzeros(TN(:,:,n)));
        spec(n,1) = TN_all(n)/(TN_all(n)+FP_all(n));
        sens(n,1) = TP_all(n)/(TP_all(n)+FN_all(n));
        ppv(n,1) = TP_all(n)/(TP_all(n)+FP_all(n));
        npv(n,1) = TN_all(n)/(TN_all(n)+FN_all(n));
        d_roc(n,1) = sqrt((1-sens(n))^2+ (1-spec(n))^2);
        d_prc(n,1) = sqrt((1-sens(n))^2+ (1-ppv(n))^2);  
        F_score(n,1) = TP_all(n)/(TP_all(n)+0.5*(FP_all(n)+FN_all(n)));
        end

% analyze the performances per subject
TP_subj = NaN(1);
FN_subj = NaN(1);
FP_subj = NaN(1);
TN_subj = NaN(1);
spec_subj = NaN(1); % specificity
sens_subj = NaN(1); % sensitivity
ppv_subj = NaN(1); % positive predictive value
npv_subj = NaN(1); % negative predicitive value
d_roc_subj = NaN(1); % distance to top left corner of ROC 
d_prc_subj = NaN(1); % distance to top left corner of PRC (precision-recall curve

 for subj = 1:size(dataBase,2)
       for n = 1:length(combination)
        TP_subj(n,subj) = sum(nonzeros(TP(subj,:,n)));
        FN_subj(n,subj) = sum(nonzeros(FN(subj,:,n)));
        FP_subj(n,subj) = sum(nonzeros(FP(subj,:,n)));
        TN_subj(n,subj) = sum(nonzeros(TN(subj,:,n)));
        spec_subj(n,subj) = TN_subj(n,subj)/(TN_subj(n,subj)+FP_subj(n,subj));
        sens_subj(n,subj) = TP_subj(n,subj)/(TP_subj(n,subj)+FN_subj(n,subj));
        ppv_subj(n,subj) = TP_subj(n,subj)/(TP_subj(n,subj)+FP_subj(n,subj));
        npv_subj(n,subj) = TN_subj(n,subj)/(TN_subj(n,subj)+FN_subj(n,subj));
        d_roc_subj(n,subj) = sqrt((1-sens_subj(n,subj))^2+ (1-spec_subj(n,subj))^2);
        d_prc_subj(n,subj) = sqrt((1-sens_subj(n,subj))^2+ (1-ppv_subj(n,subj))^2);     
       end
 end

% best combination for this run
[value, best_comb] = min(d_prc);
best_amplTh = combination(best_comb,1);
best_minSD = combination(best_comb,2);
best_sel = combination(best_comb,3);
%% more sensitive algorithm
% best combination for this run
[value2, best_comb2] = min(d_roc);
best_amplTh2 = combination(best_comb2,1);
best_minSD2 = combination(best_comb2,2);
best_sel2 = combination(best_comb2,3);
%% algorithm based on fscore
% best combination for this run
[value3, best_comb3] = min(F_score);
best_amplTh3 = combination(best_comb3,1);
best_minSD3 = combination(best_comb3,2);
best_sel3 = combination(best_comb3,3);
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
print(f1,'-dpng', 'ROC_5','-r300')

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
print(f2,'-dpng', 'PRC_5','-r300')

%% test prut (kan er later uit)
best_FP = FP_all(best_comb2,1);
best_FN = FN_all(best_comb2,1);
best_FP_subj = FP_subj(best_comb2,:);
best_FN_subj = FN_subj(best_comb2,:);
all_subj = FN_subj(best_comb,:)+FP_subj(best_comb,:)+TN_subj(best_comb,:)+TP_subj(best_comb,:);
FN_norm_subj =best_FN_subj./all_subj;
FP_norm_subj =best_FP_subj./all_subj;
%% overview of performance values of best combination for this run
matrix = zeros(6,4);

matrix(1,2:4) = spec_subj(best_comb,:);
matrix(2,2:4) =  sens_subj(best_comb,:);
matrix(3,2:4) =  ppv_subj(best_comb,:);
matrix(4,2:4) =  npv_subj(best_comb,:);
matrix(5,2:4) =  d_roc_subj(best_comb,:);
matrix(6,2:4) =  d_prc_subj(best_comb,:);

matrix(1,1) = spec(best_comb,:);
matrix(2,1) =  sens(best_comb,:);
matrix(3,1) =  ppv(best_comb,:);
matrix(4,1) =  npv(best_comb,:);
matrix(5,1) =  d_roc(best_comb,:);
matrix(6,1) =  d_prc(best_comb,:);

matrix = matrix.*100;
matrix = round(matrix);

%%
matrix2 = zeros(6,4);

matrix2(1,2:4) = spec_subj(best_comb2,:);
matrix2(2,2:4) =  sens_subj(best_comb2,:);
matrix2(3,2:4) =  ppv_subj(best_comb2,:);
matrix2(4,2:4) =  npv_subj(best_comb2,:);
matrix2(5,2:4) =  d_roc_subj(best_comb2,:);
matrix2(6,2:4) =  d_prc_subj(best_comb2,:);

matrix2(1,1) = spec(best_comb2,:);
matrix2(2,1) =  sens(best_comb2,:);
matrix2(3,1) =  ppv(best_comb2,:);
matrix2(4,1) =  npv(best_comb2,:);
matrix2(5,1) =  d_roc(best_comb2,:);
matrix2(6,1) =  d_prc(best_comb2,:);

matrix2 = matrix2.*100;
matrix2 = round(matrix2);
%%

best_spec_subj = spec_subj(best_comb,:);
best_sens_subj =  sens_subj(best_comb,:);
best_ppv_subj =  ppv_subj(best_comb,:);
best_npv_subj =  npv_subj(best_comb,:);
best_d_roc_subj =  d_roc_subj(best_comb,:);
best_d_prc_subj =  d_prc_subj(best_comb,:);

best_spec = spec(best_comb,:);
best_sens =  sens(best_comb,:);
best_ppv =  ppv(best_comb,:);
best_npv =  npv(best_comb,:);
best_d_roc =  d_roc(best_comb,:);
best_d_prc =  d_prc(best_comb,:);


%% show falsely detected cceps
% plot false positives 
cfg.amplitude_thresh = best_amplTh;
cfg.minSD = best_minSD;
cfg.sel = best_sel;
dataBase = detect_n1peak_ccep(dataBase, cfg);
CCEP_FP = NaN(1); % number of FP a CCEP by second observation
per_CCEP_FP = NaN(1);
for subj = 1:size(dataBase,2)
    for run = 1:size(dataBase(subj).metadata,2)
            detected = dataBase(subj).metadata(run).ccep.n1_peak_sample; 
            detected(~isnan(detected)) = 1;
            detected(isnan(detected)) = 0;
            scored = dataBase_visualScores(subj).metadata(run).visual_scored;
            [FP_chan,FP_stimp]= find( scored == 0 & detected == 1);
            data_FP = show_ccep(dataBase,FP_stimp,FP_chan,subj,run,cfg);
            dataBase(subj).metadata(run).ccep.n1_peak_sample_FP = data_FP;
            CCEP_FP(subj,run) = sum(data_FP(:)==1); % de CCEPs die eigenlijk wel TP zijn
            per_CCEP_FP(subj,run) = CCEP_FP(subj,run)/length(FP_chan)*100;
            clear('data_FP')
    end
end
%%
% plot false negatives
CCEP_FN = NaN(1); % number of FP a CCEP by second observation
per_CCEP_FN = NaN(1);
for subj = 1:size(dataBase,2)
    for run = 1:size(dataBase(subj).metadata,2)
            detected = dataBase(subj).metadata(run).ccep.n1_peak_sample; 
            detected(~isnan(detected)) = 1;
            detected(isnan(detected)) = 0;
            scored = dataBase_visualScores(subj).metadata(run).visual_scored;
            [FN_chan,FN_stimp]= find( scored == 1 & detected == 0); 
            data_FN = show_ccep(dataBase,FN_stimp,FN_chan,subj,run,cfg);
            dataBase(subj).metadata(run).ccep.n1_peak_sample_FN = data_FN;
            CCEP_FN(subj,run) = sum(data_FN(:)==1); % de CCEPs die echt FN zijn
            per_CCEP_FN(subj,run) = CCEP_FN(subj,run)/length(FN_chan)*100;
            clear('data_FN')
    end
end

%% peakfinder step for step
cfg.amplitude_thresh = 3.8;
cfg.minSD = 18;
cfg.sel = 14;
dataBase = detect_n1peak_ccep(dataBase, cfg);

amplitude_thresh = cfg.amplitude_thresh;
n1_peak_range = cfg.n1_peak_range;
epoch_prestim = cfg.epoch_prestim;
epoch_length = cfg.epoch_length;
minSD = cfg.minSD;
sel = cfg.sel;

subj=1;
run=1;
signal = dataBase(subj).metadata(run).cc_epoch_sorted_reref_avg;
chan=39;
stimp=2;
% create time struct 
tt = (1:epoch_length*dataBase(subj).metadata(run).ccep_header.Fs) / ...
      dataBase(subj).metadata(run).ccep_header.Fs - epoch_prestim;
               
% baseline subtraction: take median of part of the averaged signal for
% this stimulation pair before stimulation, which is the half of the
% epoch                
baseline_tt = tt>-2 & tt<-.1;
signal_median = median(signal(chan,stimp,baseline_tt),3,'omitnan');
% subtract median baseline from signal
new_signal = squeeze(signal(chan,stimp,:)) - signal_median;
peakfinder(new_signal(find(tt>0,1):find(tt>0.5,1)),sel,[],-1); % negative peaks
peakfinder(new_signal(find(tt>0,1):find(tt>0.5,1)),sel,[],1); % positive peaks
%% determine optimal settings

