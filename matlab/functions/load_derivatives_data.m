% function to read derivatives data into dataBase
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_derivatives_data(myDataPath,cfg)

input = myDataPath.input;
input_dev = myDataPath.input_dev;

dataBase = struct([]);

for nSubj = 1:size(cfg.sub_label,2)
    sub_label = ['sub-' cfg.sub_label{nSubj}];
    ses_label = cfg.ses_label{nSubj};
    task_label = cfg.task_label{nSubj};
    run_label = cfg.run_label{nSubj};

    % load electrodes.tsv and detected ccep data
    for nRun = 1:size(run_label,2)

        % load electrodes
        D = dir(fullfile(input,sub_label,ses_label,'ieeg',...
            [sub_label '_' ses_label ,'_electrodes.tsv']));

        elecsName = fullfile(D(1).folder, D(1).name);

        % remove electrodes named 'other' (so no grid, depth,strip)
        tb_electrodes = readtable(elecsName,'FileType','text','Delimiter','\t');
        idx_elec_incl = ~strcmp(tb_electrodes.group,'other');
        tb_electrodes = tb_electrodes(idx_elec_incl,:);

        % load detected and visually checked CCEP-N1s
        DD = dir(fullfile(input_dev,sub_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{nRun},'_N1sChecked.mat']));

        dataName = fullfile(DD(1).folder, DD(1).name);
        ccep = load(dataName);

        % load channels
        channelsName = fullfile(D(1).folder, DD(1).name);
        channelsName = replace(channelsName,'_N1sChecked.mat','_channels.tsv');

        tb_channels = readtable(channelsName,'FileType','text','Delimiter','\t');
        % remove electrodes not type ECOG or SEEG
        idx_ch_incl = strcmp(tb_channels.type,'ECOG')|strcmp(tb_channels.type,'SEEG');

        tb_channels = tb_channels(idx_ch_incl,:);

        dataBase(nSubj).sub_label = sub_label;
        dataBase(nSubj).ses_label = ses_label;
        dataBase(nSubj).task_label = task_label;
        dataBase(nSubj).tb_electrodes = tb_electrodes;
        dataBase(nSubj).metadata(nRun).tb_channels = tb_channels;
        dataBase(nSubj).metadata(nRun).run_label = run_label{nRun};
        dataBase(nSubj).metadata(nRun).ccep = ccep;
        fprintf('...Subject %s %s has been run...\n',sub_label,run_label{nRun})

    end
end

disp('All subjects are loaded')

