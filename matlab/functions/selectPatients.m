function cfg = selectPatients(cfg, myDataPath)

% INPUT:
% cfg.sub_label = {'RESPXXXX', 'RESPXXXX'}
% myDataPath.dataPath = string with folder location

% OUTPUT:
% cfg.sub_label = cell([1 x number of subj]), {'RESPXXXX'} {'RESPXXXX'}
% cfg.ses_label = cell([1 x number of subj]), automatically set if there is
%                   only 1 session, otherwise, 1 session can be chosen
% cfg.task_label = cell([1 x number of subj]), default SPESclin -->
%                   {'task-SPESclin'}
% cfg.run_label = cell([1 x number of subj]), each cell contains cell([1 x
%                   number of runs]) --> {'run-xxxxxx'} {'run-xxxxxx'}
%%
% for each subject
for subj = 1:size(cfg.sub_label,2)

    % find all session of this subject
    files = dir(fullfile(myDataPath.dataPath,['sub-' cfg.sub_label{subj}]));
    idx = contains({files(:).name},{'ses-'});
    sesfiles = {files(idx).name};
    string = [repmat('%s, ',1,size(sesfiles,2)-1),'%s'];

    % select one session if only one session is available
    if size(sesfiles,1) <1
        error('Pathname %s does not contain any files',fullfile(myDataPath.dataPath,['sub-' cfg.sub_label{subj}]))
    elseif size(sesfiles,1) == 1
        cfg.ses_label{subj} = sesfiles{1};
    else
        cfg.ses_label{subj} = input(sprintf(['Select one of these sessions [',string,']: \n'],sesfiles{:}),'s');
    end
end

% select task-label
task_label = cell(size(cfg.sub_label));
[task_label{:}] = deal('task-SPESclin');
cfg.task_label = task_label;

% select run: display strings to choose from
for subj = 1:size(cfg.sub_label,2)
    files = dir(fullfile(myDataPath.dataPath,['sub-' cfg.sub_label{subj}],cfg.ses_label{subj},'ieeg'));
    idx = contains({files(:).name},{'task-SPESclin'}) & contains({files(:).name},{'.eeg'});
    runfiles = {files(idx).name};
    runfiles = extractBetween(runfiles,'clin_','_ieeg');
    string = [repmat('%s, ',1,size(runfiles,2)-1),'%s'];
    
    if size(files,1) <1
        error('Pathname does not contain any ieeg-files')
    elseif size(runfiles,2) == 1
        cfg.run_label{subj} = runfiles(1);
    else
        x = input(sprintf('Select al runs of subject %s? (y/n): ',cfg.sub_label{subj}),'s'); 
        if strcmp(x,'y') 
        cfg.run_label{subj} = runfiles;
        else
        run_label = input(sprintf(['Subject %s: Select one of these files, or type more runs by separating them with a comma: ' string, ': \n'],cfg.sub_label{subj},runfiles{:}),'s');
        cfg.run_label{subj} = strsplit(run_label,{', ',','});
        end
    end
end
