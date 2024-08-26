function cfg = selectPatients(cfg, myDataPath)

% INPUT:
% cfg.sub_label = {'STREEFXX', 'STREEFXX'}
% myDataPath.input = string with folder location

% OUTPUT:
% cfg.sub_label = cell([1 x number of subj]), {'STREEFXX'} {'STREEFXX'}
% cfg.ses_label = cell([1 x number of subj]), automatically set if there is
%                   only 1 session, otherwise, 1 session can be chosen
% cfg.task_label = cell([1 x number of subj]), default SPESclin -->
%                   {'task-SPESclin'}
% cfg.run_label = cell([1 x number of subj]), each cell contains cell([1 x
%                   number of runs]) --> {'run-xxxxxx'} {'run-xxxxxx'}

%%
% for each subject

% select ses-label
ses_label = cell(size(cfg.sub_label));
[ses_label{:}] = deal('ses-1');
cfg.ses_label = ses_label;

% select task-label
task_label = cell(size(cfg.sub_label));
[task_label{:}] = deal('task-SPESclin');
cfg.task_label = task_label;

% select run
for nSubj = 1:size(cfg.sub_label,2)
    files = dir(fullfile(myDataPath.input,['sub-' cfg.sub_label{nSubj}],cfg.ses_label{nSubj},'ieeg'));
    idx = contains({files(:).name},{'task-SPESclin'}) & contains({files(:).name},{'.eeg'});
    runfiles = {files(idx).name};
    runfiles = extractBetween(runfiles,'clin_','_ieeg');

    cfg.run_label{nSubj} = runfiles;

end
end
