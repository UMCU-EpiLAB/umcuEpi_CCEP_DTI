% function to read derivatives data into dataBase
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_derivatives_data(myDataPath,cfg)

input = myDataPath.input;
input_dev = myDataPath.input_dev;

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
       

    % load electrodes.tsv and detected ccep data 
    for j=1:size(run_label,2)
        D = dir(fullfile(input,sub_label,ses_label,'ieeg',...
            [sub_label '_' ses_label '_' task_label ,'_',run_label{j}, '_ieeg.eeg']));
        
        if size(D,1) == 0
            error('%s does not exist',fullfile(input,sub_label,ses_label,'ieeg',...
                [sub_label '_' ses_label '_' task_label ,'_',run_label{j}, '_ieeg.eeg']))
        end
        
        
        % load electrodes
        D = dir(fullfile(input,sub_label,ses_label,'ieeg',...
            [sub_label '_' ses_label ,'_electrodes.tsv']));
        
        elecsName = fullfile(D(1).folder, D(1).name);
        
        % remove electrodes named 'other' (so no grid, depth,strip)
        tb_electrodes = readtable(elecsName,'FileType','text','Delimiter','\t');
        idx_elec_incl = ~strcmp(tb_electrodes.group,'other');
        tb_electrodes = tb_electrodes(idx_elec_incl,:);
         
        DD = dir(fullfile(input_dev,sub_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sChecked.mat'])); 
        if size(DD,1) == 0
            error('%s does not exist',fullfile(input_dev,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sChecked.mat']));
        end 
        dataName= fullfile(DD(1).folder, DD(1).name);
        ccep = load(dataName);

         % load channels
        channelsName = fullfile(D(1).folder, DD(1).name);
        channelsName = replace(channelsName,'_N1sChecked.mat','_channels.tsv');
        
        tb_channels = readtable(channelsName,'FileType','text','Delimiter','\t');
        % remove electrodes not type ECOG or SEEG
        idx_ch_incl = strcmp(tb_channels.type,'ECOG')|strcmp(tb_channels.type,'SEEG');
        
        tb_channels = tb_channels(idx_ch_incl,:);

        dataBase(i).sub_label = sub_label;
        dataBase(i).ses_label = ses_label;
        dataBase(i).task_label = task_label;
        dataBase(i).tb_electrodes = tb_electrodes;
        dataBase(i).ccep(j).tb_channels = tb_channels;
        dataBase(i).ccep(j).run_label = run_label{j};
        dataBase(i).ccep(j).ccep = ccep;
        fprintf('...Subject %s %s has been run...\n',sub_label,run_label{j}) 

    end
end

disp('All subjects are loaded')

