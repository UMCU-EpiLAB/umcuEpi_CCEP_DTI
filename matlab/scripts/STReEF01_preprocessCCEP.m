% STReEF01_preprocessCCEP
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 

% author: Dorien van Blooijs & Susanne Jelsma
% date: October 2021

% This script:
% - set correct paths
% - define subjects to include
% - load ECoG data
% - split into stimulation trials, 
% - re-reference the data,
% - detect CCEPs
% - visually check the detected CCEPs
% - save pre-processed data in derivatives for further analyses

%% set paths
% set umcuEpi_CCEP_DTI/matlab in your directory and run this section

clc
clear

rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

myDataPath = STReEF_setLocalDataPath(1);

% housekeeping 
clear rootPath RepoPath matlabFolder

%% patient characteristics

sub_label = ['STREEF01, STREEF02, STREEF03, STREEF04, STREEF05, STREEF06,' ...
    'STREEF07, STREEF08, STREEF09, STREEF10, STREEF11, STREEF12, STREEF13'];

cfg.sub_label = strsplit(sub_label,{', ',','});
cfg = selectPatients(cfg, myDataPath);

%% load ECoGs 

dataBase = load_ECoGdata(myDataPath,cfg);

disp('All ECoGs are loaded')

%% preprocessing CCEP in ECoG

% select epochs and average
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

dataBase = preprocess_ECoG(dataBase, cfg);

disp('All ECoGs are preprocessed')

%% rereference data
% find 10 signals in the same trial with lowest variance and not being a
% bad channel and use that as reference signal to increase the
% signal-to-noise ratio

cfg.reref = 'y';

for nSubj = 1:size(dataBase,2)
    for nRun = 1:size(dataBase(nSubj).metadata,2)
        
        dataBase = rerefCCEPdata(dataBase,nSubj,nRun,cfg);
        
        % take mean of these re-referenced signals
        dataBase(nSubj).metadata(nRun).cc_epoch_sorted_reref_avg = squeeze(mean(dataBase(nSubj).metadata(nRun).cc_epoch_sorted_reref,2,'omitnan')); 
    
        fprintf('...%s %s has been re-referenced... \n',...
            dataBase(nSubj).sub_label,dataBase(nSubj).metadata(nRun).run_label)

    end
end

disp('Data is rereferenced and all trials of one stimulus pair are averaged.')

%% detect CCEPs
% set correct parameters for detection of CCEPs in sEEG/ECoG

cfg.n1_peak_range = 100;

for nSubj = 1:size(dataBase,2)
    if all(strcmpi(dataBase(nSubj).metadata(1).tb_channels.type,'SEEG')==1)
        cfg.amplitude_thresh = 3.5;
        cfg.minSD = 16;
        cfg.sel = 0;
        dataBase(nSubj) = detect_n1peak_sEEG_ccep(dataBase(nSubj), cfg);
        fprintf('%s: CCEPs detected in sEEG \n', dataBase(nSubj).sub_label)

    elseif all(strcmpi(dataBase(nSubj).metadata(1).tb_channels.type,'ECOG')==1)
        cfg.amplitude_thresh = 2.6;
        cfg.minSD = 50;
        cfg.sel = 20;
        dataBase(nSubj) = detect_n1peak_ECoG_ccep(dataBase(nSubj), cfg);
        fprintf('%s: CCEPs detected in subdural ECoG \n', dataBase(nSubj).sub_label)
    end
end

disp('All CCEPs are detected')

%% visually check detected ccep
%check only the detected ccep's!
close all

% select subject
subs = {dataBase(:).sub_label};
string = [repmat('%s, ',1,size(subs,2)-1), '%s'];
substring = input(sprintf(['Choose subject: ',string,'\n'],subs{:}),'s');
subj = find(contains({dataBase(:).sub_label},substring));

if isempty(subj)
   error('No present subject was selected') 
end
% select run if multiple runs
if size(dataBase(subj).metadata,2)>1
    runs = {dataBase(subj).metadata(:).run_label};
    string = [repmat('%s, ',1,size(runs,2)-1), '%s'];
    substring = input(sprintf(['Choose run: ',string,'\n'],runs{:}),'s');
    run = find(contains({dataBase(subj).metadata(:).run_label},substring));
    
    if isempty(run)
        error('No present run was selected')
    end
else
    run = 1;
end

% load checked N1s if visual rating has started earlier
if exist(fullfile(myDataPath.input_dev, dataBase(subj).sub_label, ...
        [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',...
        dataBase(subj).metadata(run).task_label,'_',...
        dataBase(subj).metadata(run).run_label,'_N1sChecked.mat']),'file')

    dataBase(subj).metadata(run).ccep = load(fullfile(myDataPath.input_dev, ...
        dataBase(subj).sub_label,...
        [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',...
        dataBase(subj).metadata(run).task_label,'_',...
        dataBase(subj).metadata(run).run_label,'_N1sChecked.mat']));

end

% continue with the stimulation pair after the last saved stimulation pair
if any(strcmp(fieldnames(dataBase(subj).metadata(run)),'ccep'))
    if any(contains(fieldnames(dataBase(subj).metadata(run).ccep),'checkUntilStimp'))
        endstimp = dataBase(subj).metadata(run).ccep.checkUntilStimp;
    else
        endstimp = 0;
    end
else
    endstimp = 0;
end

% perform the visual rating
dataBase = visualRating_ccep(dataBase,myDataPath,subj,run,cfg,endstimp);

% write down whether you only checked the presence of an ER (write down
% ER), or also the location of the N1-peak (write down N1)
cfg.check_ER_N1 = input('Did you only check the presence of a CCEP, or also the location of the N1-peak? [CCEP/N1]: ','s');

% housekeeping 
clear endstimp string sub_label subs substring subj run

%% save ccep

for nSubj = 1:size(dataBase,2)
    for nRun = 1:size(dataBase(nSubj).metadata,2)
        targetFolder = fullfile(myDataPath.output, dataBase(nSubj).sub_label);
        
        % Create the folder if it doesn't exist already.
        if ~exist(targetFolder, 'dir')
            mkdir(targetFolder);
        end
        
        [~,filenameTmp,~] = fileparts(dataBase(nSubj).metadata(nRun).dataName);
        fileName = [replace(filenameTmp,'_ieeg',''),'_N1sChecked.mat'];
                
        ccep = dataBase(nSubj).metadata(nRun).ccep;
        ccep.dataName = dataBase(nSubj).metadata(nRun).dataName;
        ccep.ch = dataBase(nSubj).metadata(nRun).ch;
        ccep.cc_stimchans = dataBase(nSubj).metadata(nRun).cc_stimchans;
        ccep.cc_stimsets = dataBase(nSubj).metadata(nRun).cc_stimsets;
        ccep.epoch_length = cfg.epoch_length;
        ccep.epoch_prestim = cfg.epoch_prestim;
        ccep.reref = cfg.reref;
        
        save(fullfile(targetFolder,fileName), '-struct','ccep');
        fprintf('Saved ccep-struct in %s \n',fullfile(targetFolder,fileName))

    end
end
