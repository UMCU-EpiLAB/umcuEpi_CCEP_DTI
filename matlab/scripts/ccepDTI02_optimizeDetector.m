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

clear sub_label sub_label1 sub_label2 files1 files2 i x

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
                data_sorted_avg = dataBase(subj).metadata_runs(run).cc_epoch_sorted_avg; 
                data_sorted_reref_avg = dataBase(subj).metadata_runs(run).cc_epoch_sorted_reref_avg;
                stimsets = dataBase(subj).metadata_runs(run).cc_stimsets;
    
    
                dataBase(subj).metadata.ccep_VS1.checked = [dataBase(subj).metadata.ccep_VS1.checked scored];
                dataBase(subj).metadata.ccep_VS2.checked = [dataBase(subj).metadata.ccep_VS2.checked scored2];
                dataBase(subj).metadata.cc_epoch_sorted_avg = [dataBase(subj).metadata.cc_epoch_sorted_avg data_sorted_avg];
                dataBase(subj).metadata.cc_epoch_sorted_reref_avg = [dataBase(subj).metadata.cc_epoch_sorted_reref_avg data_sorted_reref_avg];
                dataBase(subj).metadata.cc_stimsets = [dataBase(subj).metadata.cc_stimsets;stimsets];
                %voeg hier nog alle andere data die je gebruikt die per run
                %verschilt toe!!!
            end
    end
end
clear data data_sorted_reref_avg data_sorted_avg run stimsets
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
clear subjs subj scored2 scored
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
%% optimize detector
% for cECoG data, these are the best parameters:
cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;
cfg.minSD = 50;
cfg.sel = 20;

% range of parameters tested:
minSD_range = 11:1:13 ;
sel_range = 10:1:12 ;
amplitude_tresh_range = 3.2:0.1:3.3 ;

n=1;
% pre-allocation
subjs = size(dataBase,2);
combs = size(minSD_range,2)*size(sel_range,2)*size(amplitude_tresh_range,2);
TP = NaN(combs,subjs); % true positives
FN = NaN(combs,subjs); % false negatives
FP = NaN(combs,subjs); % false positives
TN = NaN(combs,subjs); % true negatives
combination = NaN(combs,3); % combination of parameters. dit geeft een error

for amplTh = amplitude_tresh_range 
for sd = minSD_range 
for sl = sel_range
    cfg.amplitude_thresh = amplTh;
    cfg.minSD = sd;
    cfg.sel = sl;
    dataBase = detect_n1peak_ccep(dataBase, cfg); 
    combination(n,:) = [amplTh sd sl]; % combination(n= number combination, parameter) parameter 1= amplitude_thresh; 2=minSD; 3=sel
    for subj = 1:size(dataBase,2)
        detected = dataBase(subj).metadata.ccep.n1_peak_sample; 
        detected(~isnan(detected)) = 1;
        detected(isnan(detected)) = 0;
        scored = dataBase(subj).metadata.visual_scored;

        TP(n,subj) = numel(find( scored == 1 & detected == 1)); 
        FN(n,subj) = numel(find( scored == 1 & detected == 0));
        FP(n,subj) = numel(find( scored == 0 & detected == 1));
        TN(n,subj) = numel(find( scored == 0 & detected == 0));
    end
    n=n+1;
    fprintf('...Combination: Amplitude Threshold=%g, min SD of baseline=%g, Sel=%g, has been run...\n',...
        amplTh,sd,sl)
end
end
end
clear n subjs combs sel_range minSD_range amplitude_tresh_range detected scored sd sl amplTh

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
    else
        error('Number of visual scored CCEPs is not the same as sum of all detected');
    end
end
clear subjs subj n 
%% best combination for this range of parameters
[value, best_comb] = min(d_prc);
best_amplTh = combination(best_comb,1);
best_minSD = combination(best_comb,2);
best_sel = combination(best_comb,3);

%% more sensitive algorithm
% best combination for this range of parameters
[value2, best_comb2] = min(d_roc);
best_amplTh2 = combination(best_comb2,1);
best_minSD2 = combination(best_comb2,2);
best_sel2 = combination(best_comb2,3);
%% algorithm based on fscore
% best combination for this range of parameters
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
%saveas(f1,'ROC_5_d_roc.m')

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
% best_FP = FP_all(best_comb2,1);
% best_FN = FN_all(best_comb2,1);
% best_FP_subj = FP_subj(best_comb2,:);
% best_FN_subj = FN_subj(best_comb2,:);
% all_subj = FN_subj(best_comb,:)+FP_subj(best_comb,:)+TN_subj(best_comb,:)+TP_subj(best_comb,:);
% FN_norm_subj =best_FN_subj./all_subj;
% FP_norm_subj =best_FP_subj./all_subj;
%% overview of performance values of best combination for this range of parameters
%je zet dit nu in een matrix, wat prima is, maar je zou het ook in een tabel kunnen zetten 
% en dit dan kunnen weergeven in je command window. 
% Dan heb je direct je resultaten duidelijk in beeld. 
% Ik weet in de matrix ook niet welke volgorde van variabelen (sens, spec etc) je gebruikt, 
% in een tabel kan dit in 1x overzichtelijk zijn.
pf_value = { 'Specificity';'Sensitivity';'Positive predicitve value';'Negative predicitve value';'Distance to corner ROC (receiver operating characteristic) curve';'Distance to corner precision-recall curve (PRC)'};
pf_overal = round([spec(best_comb,:); sens(best_comb,:); ppv(best_comb,:);npv(best_comb,:);d_roc(best_comb,:);d_prc(best_comb,:)].*100);
pf_subj = round([spec_subj(best_comb,:); sens_subj(best_comb,:);ppv_subj(best_comb,:);npv_subj(best_comb,:);d_roc_subj(best_comb,:); d_prc_subj(best_comb,:)].*100);
pf_table = table(pf_value,pf_overal,pf_subj(:,1),pf_subj(:,2),pf_subj(:,3),'VariableNames',{'Performance','Overall','Subject 1','Subject 2','Subject 3'});

%% show falsely detected cceps
% % plot false positives 
% cfg.amplitude_thresh = best_amplTh;
% cfg.minSD = best_minSD;
% cfg.sel = best_sel;
% dataBase = detect_n1peak_ccep(dataBase, cfg);
% CCEP_FP = NaN(1); % number of FP a CCEP by second observation
% per_CCEP_FP = NaN(1);
% for subj = 1:size(dataBase,2)
%     for run = 1:size(dataBase(subj).metadata,2)
%             detected = dataBase(subj).metadata(run).ccep.n1_peak_sample; 
%             detected(~isnan(detected)) = 1;
%             detected(isnan(detected)) = 0;
%             scored = dataBase_visualScores(subj).metadata(run).visual_scored;
%             [FP_chan,FP_stimp]= find( scored == 0 & detected == 1);
%             data_FP = show_ccep(dataBase,FP_stimp,FP_chan,subj,run,cfg);
%             dataBase(subj).metadata(run).ccep.n1_peak_sample_FP = data_FP;
%             CCEP_FP(subj,run) = sum(data_FP(:)==1); % de CCEPs die eigenlijk wel TP zijn
%             per_CCEP_FP(subj,run) = CCEP_FP(subj,run)/length(FP_chan)*100;
%             clear('data_FP')
%     end
% end
% %%
% % plot false negatives
% CCEP_FN = NaN(1); % number of FP a CCEP by second observation
% per_CCEP_FN = NaN(1);
% for subj = 1:size(dataBase,2)
%     for run = 1:size(dataBase(subj).metadata,2)
%             detected = dataBase(subj).metadata(run).ccep.n1_peak_sample; 
%             detected(~isnan(detected)) = 1;
%             detected(isnan(detected)) = 0;
%             scored = dataBase_visualScores(subj).metadata(run).visual_scored;
%             [FN_chan,FN_stimp]= find( scored == 1 & detected == 0); 
%             data_FN = show_ccep(dataBase,FN_stimp,FN_chan,subj,run,cfg);
%             dataBase(subj).metadata(run).ccep.n1_peak_sample_FN = data_FN;
%             CCEP_FN(subj,run) = sum(data_FN(:)==1); % de CCEPs die echt FN zijn
%             per_CCEP_FN(subj,run) = CCEP_FN(subj,run)/length(FN_chan)*100;
%             clear('data_FN')
%     end
% end
% 
%% peakfinder step for step
% cfg.amplitude_thresh = 3.8;
% cfg.minSD = 18;
% cfg.sel = 14;
% dataBase = detect_n1peak_ccep(dataBase, cfg);
% 
% amplitude_thresh = cfg.amplitude_thresh;
% n1_peak_range = cfg.n1_peak_range;
% epoch_prestim = cfg.epoch_prestim;
% epoch_length = cfg.epoch_length;
% minSD = cfg.minSD;
% sel = cfg.sel;
% 
% subj=1;
% run=1;
% signal = dataBase(subj).metadata(run).cc_epoch_sorted_reref_avg;
% chan=39;
% stimp=2;
% % create time struct 
% tt = (1:epoch_length*dataBase(subj).metadata(run).ccep_header.Fs) / ...
%       dataBase(subj).metadata(run).ccep_header.Fs - epoch_prestim;
%                
% % baseline subtraction: take median of part of the averaged signal for
% % this stimulation pair before stimulation, which is the half of the
% % epoch                
% baseline_tt = tt>-2 & tt<-.1;
% signal_median = median(signal(chan,stimp,baseline_tt),3,'omitnan');
% % subtract median baseline from signal
% new_signal = squeeze(signal(chan,stimp,:)) - signal_median;
% peakfinder(new_signal(find(tt>0,1):find(tt>0.5,1)),sel,[],-1); % negative peaks
% peakfinder(new_signal(find(tt>0,1):find(tt>0.5,1)),sel,[],1); % positive peaks
%% determine optimal settings

