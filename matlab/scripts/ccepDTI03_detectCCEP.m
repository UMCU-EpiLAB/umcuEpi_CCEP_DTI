% ccepDTI03_detectCCEP
% detect CCEPs in either cECoG or sEEG data

% author: Dorien van Blooijs & Susanne Jelsma
% date: October 2021

% load ECoG data, split into stimulation trials, re-reference the data,
% detect CCEPs and visually check the detected CCEPs

%% set paths
% set umcuEpi_CCEP_DTI/matlab in your directory and run this section

clc
clear
cfg.folderinput = 'chronic_ECoG'; % from which folder would you like to load ECoGs?
myDataPath = setLocalDataPath(cfg);

%% patient characteristics
sub_label = input('Patient number (RESPXXXX) (select multiple patients by separating each with a comma): ','s');

cfg.sub_label = strsplit(sub_label,{', ',','});

cfg = selectPatients(cfg, myDataPath);

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

%% detect CCEPs

% set correct parameters for detection of CCEPs in sEEG/ECoG

% TODO: in the function detect_n1peak_ccep, all subjects and runs in the
% dataBase are run, but the  next few lines for setting parameters are only
% for one subject (when you would load several subjects...), because you
% could load both sEEG and ECoG patients in one database, and then
% different parameters might be necessary.
if all(strcmpi(dataBase(subj).metadata(1).tb_channels.type,'SEEG')==1)
    % to be determined in the previous step

elseif all(strcmpi(dataBase(subj).metadata(1).tb_channels.type,'ECoG')==1) % CHECK: im not sure that it is ECoG in type in ECoG data...
    cfg.amplitude_thresh = 2.6;
    cfg.n1_peak_range = 100;
end

dataBase = detect_n1peak_ECoG_ccep(dataBase, cfg);

disp('All CCEPs are detected')


%% visually check detected ccep
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
if exist(fullfile(myDataPath.CCEPpath, dataBase(subj).sub_label, ...
        dataBase(subj).ses_label, dataBase(subj).metadata(run).run_label,...
        [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',...
        dataBase(subj).metadata(run).task_label,'_',...
        dataBase(subj).metadata(run).run_label,'_N1sChecked.mat']),'file')
    
   dataBase(subj).metadata(run).ccep = load(fullfile(myDataPath.CCEPpath, ...
       dataBase(subj).sub_label,dataBase(subj).ses_label, dataBase(subj).metadata(run).run_label,...
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

% visual rating
dataBase = visualRating_ccep(dataBase,myDataPath,subj,run,cfg,endstimp);

% write down whether you only checked the presence of an ER (write down
% ER), or also the location of the N1-peak (write down N1)
cfg.check_ER_N1 = input('Did you only check the presence of a CCEP, or also the location of the N1-peak? [CCEP/N1]: ','s');

%% save ccep

for subj = 1:size(dataBase,2)
    for run = 1:size(dataBase(subj).metadata,2)
        targetFolder = fullfile(myDataPath.CCEPpath, dataBase(subj).sub_label,dataBase(subj).ses_label,dataBase(subj).metadata(run).run_label);
        
        % Create the folder if it doesn't exist already.
        if ~exist(targetFolder, 'dir')
            mkdir(targetFolder);
        end
        
        start_filename = strfind(dataBase(subj).metadata(run).dataName,'/');
        stop_filename = strfind(dataBase(subj).metadata(run).dataName,'_ieeg');
        
        fileName = [dataBase(subj).metadata(run).dataName(start_filename(end)+1:stop_filename-1),'_N1sChecked.mat'];
                
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
