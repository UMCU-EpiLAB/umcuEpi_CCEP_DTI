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

%% optimize detector

% for cECoG data, these are the best parameters:
cfg.amplitude_thresh = 2.6;
cfg.n1_peak_range = 100;

% let's find out what are the best parameters for sEEG data
% net to the amplitude_threshold, you might want to optimize pre_stim_sd
% (line121), and sel (line 145).
% and now, only the negative peaks are taken into account, but I think,
% since sEEG is stimulated differently, that you might want to take both
% negative and positive peaks into account (in peakfinder, line 145). 

n=1;
% pre-allocation
TP = NaN(1); % true positives
FN = NaN(1); % false negatives
FP = NaN(1); % false positives
TN = NaN(1); % true negatives

for amplTh = 0.5:0.5:5 % vary amplitude threshold
    cfg.amplitude_thresh = amplTh;
    
    dataBase = detect_n1peak_ccep(dataBase, cfg);

    for subj = 1:size(dataBase,2)
        for run = 1:size(dataBase(subj).metadata,2)
            detected = dataBase(subj).metadata(run).ccep.n1_peak_sample;
            detected(~isnan(detected)) = 1;

            TP(subj,run,n) = numel(find(dataBase(subj).metadata(run).ccep.checked == 1 & ...
                detected == 1));
            FN(subj,run,n) = numel(find(dataBase(subj).metadata(run).ccep.checked == 1 & ...
                detected == 0));
            FP(subj,run,n) = numel(find(dataBase(subj).metadata(run).ccep.checked == 0 & ...
                detected == 1));
            TN(subj,run,n) = numel(find(dataBase(subj).metadata(run).ccep.checked == 0 & ...
                detected == 0));
        end
    end
    n=n+1;
    
end


%% calculate and visualize performance for each setting



%% determine optimal settings



