% function to read derivatives data into dataBase
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_detected_data(dataBase,myDataPath,cfg)

input_dev = myDataPath.input_dev;

for i=1:size(cfg.sub_label,2)
    if isfield(cfg,'sub_label')
        sub_label = ['sub-' cfg.sub_label{i}];
    else
        error('No sub-label specified');
    end
    if isfield(cfg,'ses_label')
        ses_label = cfg.ses_label{i};
    else
        error('No ses-label specified');
    end
    if isfield(cfg,'task_label')
        task_label = cfg.task_label{i};
    else
        error('No task-label specified');
    end
    if isfield(cfg,'run_label')
        run_label = cell(1);
        if size(cfg.run_label{i},2) == 1
            run_label{1} = cfg.run_label{i}{1};
        elseif size(cfg.run_label{i},2) > 1
            for j = 1:size(cfg.run_label{i},2)
                run_label{j} = cfg.run_label{i}{j}; 
            end
        end
    else
        error('No run-label specified');
    end
       

    % load detected ccep data 
    for j=1:size(run_label,2)        
        DD = dir(fullfile(input_dev,sub_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sDetected.mat'])); 
        if size(DD,1) == 0
            error('%s does not exist',fullfile(input_dev,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sDetected.mat']));
        end 
        dataName= fullfile(DD(1).folder, DD(1).name);
        ccep = load(dataName);

        dataBase(i).ccep(j).run_label = run_label{j};
        dataBase(i).ccep(j).ccep = ccep;
        fprintf('...Subject %s %s has been run...\n',sub_label,run_label{j})
    end
end
disp('All subjects are loaded')

end
