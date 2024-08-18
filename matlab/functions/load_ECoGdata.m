% function to read data into dataBase
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_ECoGdata(myDataPath,cfg)

dataPath = myDataPath.input;
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
        D = dir(fullfile(dataPath,sub_label,ses_label,'ieeg',...
            [sub_label '_' ses_label '_' task_label ,'_',run_label{j}, '_ieeg.eeg']));
        
        if size(D,1) == 0
            error('%s does not exist',fullfile(dataPath,sub_label,ses_label,'ieeg',...
                [sub_label '_' ses_label '_' task_label ,'_',run_label{j}, '_ieeg.eeg']))
        end
        
        
        dataName = fullfile(D(1).folder, D(1).name);
        
        ccep_data = ft_read_data(dataName,'dataformat','brainvision_eeg');
        ccep_header = ft_read_header(dataName);
        
        % load events        
        eventsName = replace(dataName,'_ieeg.eeg','_events.tsv');
        
        tb_events = readtable(eventsName,'FileType','text','Delimiter','\t');
        
        % load electrodes
        D = dir(fullfile(dataPath,sub_label,ses_label,'ieeg',...
            [sub_label '_' ses_label ,'_electrodes.tsv']));
        
        elecsName = fullfile(D(1).folder, D(1).name);
        
        % remove electrodes named 'other' (so no grid, depth,strip)
        tb_electrodes = readtable(elecsName,'FileType','text','Delimiter','\t');
        idx_elec_incl = ~strcmp(tb_electrodes.group,'other');
        tb_electrodes = tb_electrodes(idx_elec_incl,:);
        
        % load channels
        channelsName = replace(dataName,'_ieeg.eeg','_channels.tsv');
        
        tb_channels = readtable(channelsName,'FileType','text','Delimiter','\t');
        % remove electrodes not type ECOG or SEEG
        idx_ch_incl = strcmp(tb_channels.type,'ECOG')|strcmp(tb_channels.type,'SEEG');
        
        tb_channels = tb_channels(idx_ch_incl,:);
        ch_incl = tb_channels.name;
        
        data = ccep_data(idx_ch_incl,:);
        
        dataBase(i).sub_label = sub_label;
        dataBase(i).ses_label = ses_label;
        dataBase(i).metadata(j).task_label = task_label;
        dataBase(i).metadata(j).run_label = run_label{j};
        dataBase(i).metadata(j).dataName = dataName;
        dataBase(i).metadata(j).ccep_header = ccep_header;
        dataBase(i).metadata(j).tb_events = tb_events;
        dataBase(i).metadata(j).tb_channels = tb_channels;
        dataBase(i).tb_electrodes = tb_electrodes;
        dataBase(i).metadata(j).ch = ch_incl;
        dataBase(i).metadata(j).data = data;
        fprintf('...Subject %s %s has been run...\n',sub_label,run_label{j})
    end
end

disp('All subjects are loaded')
