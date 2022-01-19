% function to read scored data into dataBase
% author: Susanne Jelsma 
% date: November 2021


function dataBase = load_visualScores(myDataPath,cfg)

dataPath = myDataPath.CCEPpath; 
if isfield(myDataPath,'CCEPpath2')
dataPath2 = myDataPath.CCEPpath2;
else
error('myDataPath.CCEPpath2 does not exist. Make sure you add this folder for a second observer in personalDataPath.m');
end

dataBase = struct([]);

if isfield(cfg,'sub_label')==0
        error('No sub-label specified');
end

for i=1:size(cfg.sub_label,2)
    sub_label = ['sub-' cfg.sub_label{i}];
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
        D2 = dir(fullfile(dataPath2,sub_label,ses_label,run_label{j},...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sChecked.mat']));
      
        if size(D,1) == 0
            error('%s does not exist',fullfile(dataPath,sub_label,ses_label,run_label{j},...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sChecked.mat']));
        end  

        if size(D2,1) == 0
            error('%s does not exist',fullfile(dataPath2,sub_label,ses_label,run_label{j},...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sChecked.mat']));
        end  
        
        dataName = fullfile(D(1).folder, D(1).name);
        dataName2 = fullfile(D2(1).folder, D2(1).name);
        
        data= load(dataName);
        data2= load(dataName2);   
        
        dataBase(i).sub_label = sub_label;
        dataBase(i).ses_label = ses_label;
        dataBase(i).metadata(j).task_label = task_label;
        dataBase(i).metadata(j).run_label = run_label{j};
        dataBase(i).metadata(j).dataName = dataName;
        dataBase(i).metadata(j).ccep = data;
        dataBase(i).metadata(j).ccep2 = data2;
        fprintf('...Subject %s %s has been run...\n',sub_label,run_label{j})
    end
end
fprintf('All subjects from %s are loaded\n', dataPath)
fprintf('All subjects from %s are loaded\n', dataPath2)

%je zou hier nog aan toe kunnen voegen uit welke file je deze info hebt gehaald, z
% odat duidelijk is in je command window wie observer 1 en wie observer 2 is.
end

