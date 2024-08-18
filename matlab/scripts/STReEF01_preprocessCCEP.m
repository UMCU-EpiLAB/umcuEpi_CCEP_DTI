% STReEF01_preprocessCCEP
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 
% detect CCEPs in either cECoG or sEEG data

% author: Dorien van Blooijs & Susanne Jelsma
% date: October 2021

% load ECoG data, split into stimulation trials, re-reference the data,
% detect CCEPs

%% set paths
% set umcuEpi_CCEP_DTI/matlab in your directory and run this section

clc
clear

rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

cfg.folderinput = 'shareData_STReEF'; % from which folder would you like to load ECoGs?
myDataPath = setLocalDataPath(cfg);

% housekeeping 
clear rootPath RepoPath matlabFolder

%% patient characteristics
sub_label = input('Patient number (STREEFXX) (select multiple patients by separating each with a comma): ','s');

cfg.sub_label = strsplit(sub_label,{', ',','});

cfg = selectPatients(cfg, myDataPath);

%% load ECoGs 

dataBase = load_ECoGdata(myDataPath,cfg);

disp('All ECoGs are loaded')

%% preprocessing CCEP in ECoG

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

%% detect CCEPs
% set correct parameters for detection of CCEPs in sEEG/ECoG

for subj = 1:size(dataBase,2)
if all(strcmpi(dataBase(subj).metadata(1).tb_channels.type,'SEEG')==1)
    cfg.amplitude_thresh = 3.5; 
    cfg.n1_peak_range = 100;
    cfg.minSD = 16;
    cfg.sel = 0;
    dataBase(subj) = detect_n1peak_sEEG_ccep(dataBase(subj), cfg);
    disp('seeg')
elseif all(strcmpi(dataBase(subj).metadata(1).tb_channels.type,'ECOG')==1) 
    cfg.amplitude_thresh = 2.6; 
    cfg.n1_peak_range = 100;
    cfg.minSD = 50;
    cfg.sel = 20;
    dataBase(subj) = detect_n1peak_ECoG_ccep(dataBase(subj), cfg);
    disp('grid')
end
end

disp('All CCEPs are detected')

%% save ccep

for subj = 1:size(dataBase,2)
    for run = 1:size(dataBase(subj).metadata,2)
        targetFolder = fullfile(myDataPath.output, dataBase(subj).sub_label);
        
        % Create the folder if it doesn't exist already.
        if ~exist(targetFolder, 'dir')
            mkdir(targetFolder);
        end
        
        start_filename = strfind(dataBase(subj).metadata(run).dataName,'/');
        stop_filename = strfind(dataBase(subj).metadata(run).dataName,'_ieeg');
        
        fileName = [dataBase(subj).metadata(run).dataName(start_filename(end)+1:stop_filename-1),'_N1sDetected.mat'];
                
        ccep = dataBase(subj).metadata(run).ccep;
        ccep.dataName = dataBase(subj).metadata(run).dataName;
        ccep.ch = dataBase(subj).metadata(run).ch;
        ccep.cc_stimchans = dataBase(subj).metadata(run).cc_stimchans;
        ccep.cc_stimsets = dataBase(subj).metadata(run).cc_stimsets;
        ccep.epoch_length = cfg.epoch_length;
        ccep.epoch_prestim = cfg.epoch_prestim;
        
        if any(contains(fieldnames(cfg),'dir'))
            ccep.dir = cfg.dir;
            ccep.amp = cfg.amp;
        end
        
        if any(contains(fieldnames(cfg),'reref'))
            ccep.reref = cfg.reref;
        end
        
        if any(contains(fieldnames(cfg),'check_ER_N1'))
           ccep.check_ER_N1 = cfg.check_ER_N1; 
        end
        
        if any(contains(fieldnames(cfg),'n1Detected'))
            ccep.n1Detected = cfg.n1Detected;
        end

        save(fullfile(targetFolder,fileName), '-struct','ccep');
        fprintf('Saved ccep-struct in %s \n',fullfile(targetFolder,fileName))

    end
end
