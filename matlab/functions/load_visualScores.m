% function to read data into dataBase
% author: Dorien van Blooijs
% date: June 2019
% Hier maak ik een functie de de gescoorde datasets inlaad voor beide
% scorers

function dataBase = load_visualScores(myDataPath,cfg, scorer)

if scorer==1
dataPath = myDataPath.CCEPpath; % CCEPpath2 is other scorer maak dit nog eleganter
else 
dataPath = myDataPath.CCEPpath2;
end

dataBase = struct([]);

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
    
    for j=1:size(run_label,2)
        D = dir(fullfile(dataPath,sub_label,ses_label,run_label{j},...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sChecked.mat']));
      
        if size(D,1) == 0
            error('%s does not exist',fullfile(dataPath,sub_label,ses_label,run_label{j},...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sChecked.mat']));
        end  
        
        dataName = fullfile(D(1).folder, D(1).name);
        
        data= load(dataName);   
        
        dataBase(i).sub_label = sub_label;
        dataBase(i).ses_label = ses_label;
        dataBase(i).metadata(j).task_label = task_label;
        dataBase(i).metadata(j).run_label = run_label{j};
        dataBase(i).metadata(j).dataName = dataName;
        dataBase(i).metadata(j).ccep = data;
        fprintf('...Subject %s %s has been run...\n',sub_label,run_label{j})
    end
end

disp('All subjects are loaded')
