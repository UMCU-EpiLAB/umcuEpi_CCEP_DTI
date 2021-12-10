% ccepDTI02_optimizeDetector

% author: Dorien van Blooijs & Susanne Jelsma
% date: October 2021

% make sure you first have some visually scored CCEPs (ccepDTI01_scoreCCEP)
% before you continue with this script

% load ECoG data, split into stimulation trials, re-reference the data,
% 


%% set paths
% set umcuEpi_CCEP_DTI/matlab in your directory and run this section

clc
clear
cfg.folderinput = 'chronic_ECoG'; % from which folder would you like to load ECoGs?
myDataPath = setLocalDataPath(cfg);

%% patient characteristics

% I think that we should first mention which patients have visually scored
% CCEPs in sEEG data, and load these specific patients

sub_label = input('Patient number (RESPXXXX) (select multiple patients by separating each with a comma): ','s'); %(select multiple patients by separating each with a comma)

cfg.sub_label = strsplit(sub_label,{', ',','});

cfg = selectPatients(cfg, myDataPath);
%% load visual scoring results
dataBase_visualScores = load_visualScores(myDataPath,cfg,1);
dataBase_visualScores2 = load_visualScores(myDataPath,cfg,2);
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
            scored2 = dataBase_visualScores2(subj).metadata(run).ccep.checked;
            [ccep_row, ccep_col] = find( scored == 1 & scored2 == 1);
            agreed_cceps = zeros(size(scored));
            for indx=1:length(ccep_row)
            agreed_cceps(ccep_row(indx),ccep_col(indx)) = 1;
            end
            dataBase_visualScores(subj).metadata(run).visual_scored = agreed_cceps;
        end
    end
clear('agreed_cceps');
clear('ccep_row');
clear('ccep_col');
clear('indx');
clear('scored');
clear('scored2');
clear('dataBase_visualScores2');
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

for sd = 0:10:100 %10:1:50
for sl = 10:10:200 %10:1:20
for amplTh = 0.5:0.5:15 %1.5:0.1:2.5 
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

        for n = 1:length(combination)
        TP_all(n,1) = (sum(nonzeros(TP(:,:,n))));
        FN_all(n,1) = sum(nonzeros(FN(:,:,n)));
        FP_all(n,1) = sum(nonzeros(FP(:,:,n)));
        TN_all(n,1) = sum(nonzeros(TN(:,:,n)));
        spec(n,1) = TN_all(n)/(TN_all(n)+FP_all(n));
        sens(n,1) = TP_all(n)/(TP_all(n)+FN_all(n));
        ppv(n,1) = TP_all(n)/(TP_all(n)+FP_all(n));
        npv(n,1) = TN_all(n)/(TN_all(n)+FN_all(n));
        d_roc(n,1) = sqrt((1-sens(n))^2+ (1-spec(n))^2);
        d_prc(n,1) = sqrt((1-sens(n))^2+ (1-ppv(n))^2);       
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
%% show falsely detected cceps
% plot false positives 
cfg.amplitude_thresh = best_amplTh;
cfg.minSD = best_minSD;
cfg.sel = best_sel;
dataBase = detect_n1peak_ccep(dataBase, cfg);
per_CCEP_FP = NaN(1); % percentage of the FP wich is by second observation TP
for subj = 1:size(dataBase,2)
    for run = 1:size(dataBase(subj).metadata,2)
            detected = dataBase(subj).metadata(run).ccep.n1_peak_sample; 
            detected(~isnan(detected)) = 1;
            detected(isnan(detected)) = 0;
            scored = dataBase_visualScores(subj).metadata(run).visual_scored;
            [FP_chan,FP_stimp]= find( scored == 0 & detected == 1);
            data_FP = show_ccep(dataBase,FP_stimp,FP_chan,subj,run,cfg);
            dataBase(subj).metadata(run).ccep.n1_peak_sample_FP = data_FP;
            per_CCEP_FP(subj,run) = sum(data_FP(:)==1)/length(FP_chan)*100;
            clear('data_FP')
    end
end
%%
% plot false negatives
per_CCEP_FN = NaN(1); % percentage of the FN wich is by second observation TN
for subj = 1:size(dataBase,2)
    for run = 1:size(dataBase(subj).metadata,2)
            detected = dataBase(subj).metadata(run).ccep.n1_peak_sample; 
            detected(~isnan(detected)) = 1;
            detected(isnan(detected)) = 0;
            scored = dataBase_visualScores(subj).metadata(run).visual_scored;
            [FN_chan,FN_stimp]= find( scored == 1 & detected == 0); 
            data_FN = show_ccep(dataBase,FN_stimp,FN_chan,subj,run,cfg);
            dataBase(subj).metadata(run).ccep.n1_peak_sample_FN = data_FN;
            per_CCEP_FN(subj,run) = (1-sum(data_FN(:)==1)/length(FN_chan))*100;
            clear('data_FN')
    end
end
%% determine optimal settings



