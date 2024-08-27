% function to read data into dataBase
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_ECoGdata(myDataPath,cfg)

dataPath = myDataPath.input;
dataBase = struct([]);

for nSubj = 1:size(cfg.sub_label,2)
    sub_label = ['sub-' cfg.sub_label{nSubj}];
    ses_label = cfg.ses_label{nSubj};
    task_label = cfg.task_label{nSubj};
    
    run_label = cell(1);
    for nRun = 1:size(cfg.run_label{nSubj},2)
        run_label{nRun} = cfg.run_label{nSubj}{nRun};
    end
    
    for nRun = 1:size(run_label,2)
        D = dir(fullfile(dataPath,sub_label,ses_label,'ieeg',...
            [sub_label '_' ses_label '_' task_label ,'_',run_label{nRun}, '_ieeg.eeg']));
        
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
        
        dataBase(nSubj).sub_label = sub_label;
        dataBase(nSubj).ses_label = ses_label;
        dataBase(nSubj).metadata(nRun).task_label = task_label;
        dataBase(nSubj).metadata(nRun).run_label = run_label{nRun};
        dataBase(nSubj).metadata(nRun).dataName = dataName;
        dataBase(nSubj).metadata(nRun).ccep_header = ccep_header;
        dataBase(nSubj).metadata(nRun).tb_events = tb_events;
        dataBase(nSubj).metadata(nRun).tb_channels = tb_channels;
        dataBase(nSubj).tb_electrodes = tb_electrodes;
        dataBase(nSubj).metadata(nRun).ch = ch_incl;
        dataBase(nSubj).metadata(nRun).data = data;
        fprintf('...Subject %s %s has been run...\n',sub_label,run_label{nRun})
    end
end

disp('All subjects are loaded')
