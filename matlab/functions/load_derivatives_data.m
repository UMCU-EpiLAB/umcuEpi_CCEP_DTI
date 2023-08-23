% function to read data into dataBase
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_derivatives_data(myDataPath,cfg)

dataPath = myDataPath.dataPath;
dataPath2 = myDataPath.DWIMATLABpath;
dataPath3 = myDataPath.CCEP;

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
    run_label_dwi = cfg.run_label{i}{1}; % er is maar een dwi per subj, zo is het opgeslagen beetje onhandig
    
    % load dwiInfo which contains information about the included electrode contact coordinates and the electrode contact areas
    DWI = dir(fullfile(dataPath2,sub_label,ses_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label_dwi,'_dwiInfo.mat'])); % made by running STReEF02_coreg_roidef_matlab

    if size(DWI,1) == 0
            error('%s does not exist',fullfile(dataPath2,sub_label,ses_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label_dwi,'_dwiInfo.mat']));
    end 
    dataNameDWI = fullfile(DWI(1).folder, DWI(1).name);
    dwi = load(dataNameDWI);

     % load networkInfo which contains information about the included electrode contact coordinates and the electrode contact areas
    NETWORK = dir(fullfile(dataPath2,sub_label,ses_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label_dwi,'_networkInfo.mat'])); % made by running STReEF05_compare_networks

    if size(NETWORK,1) == 0
            warning('%s does not exist (no problem if you are executing STReEF05)',fullfile(dataPath2,sub_label,ses_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label_dwi,'_networkInfo.mat'])); % load only if present because function also used for start of STReEF05_compare_networks 
    else
        dataNameNETWORK = fullfile(NETWORK(1).folder, NETWORK(1).name);
        network = load(dataNameNETWORK);
        dataBase(i).network = network;
    end 
    

    % load electrodes.tsv and detected and visual scored ccep data (made by running umcuEpi_CCEP_DTI (main)/matlab/scripts/ccepDTI03_detectCCEP.m) 
    for j=1:size(run_label,2)
        D = dir(fullfile(dataPath,sub_label,ses_label,'ieeg',...
            [sub_label '_' ses_label '_' task_label ,'_',run_label{j}, '_ieeg.eeg']));
        
        if size(D,1) == 0
            error('%s does not exist',fullfile(dataPath,sub_label,ses_label,'ieeg',...
                [sub_label '_' ses_label '_' task_label ,'_',run_label{j}, '_ieeg.eeg']))
        end
        
        
        dataName = fullfile(D(1).folder, D(1).name);

        % load electrodes
        D = dir(fullfile(dataPath,sub_label,ses_label,'ieeg',...
            [sub_label '_' ses_label ,'_electrodes.tsv']));
        
        elecsName = fullfile(D(1).folder, D(1).name);
        
        % remove electrodes named 'other' (so no grid, depth,strip)
        tb_electrodes = readtable(elecsName,'FileType','text','Delimiter','\t');
        idx_elec_incl = ~strcmp(tb_electrodes.group,'other');
        tb_electrodes = tb_electrodes(idx_elec_incl,:);
         
        DRE = dir(fullfile(dataPath3,sub_label,ses_label,run_label{j},...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sChecked.mat'])); % checked data
        if size(DRE,1) == 0
            error('%s does not exist',fullfile(dataPath3,sub_label,ses_label,run_label{j},...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label{j},'_N1sChecked.mat']));
        end 
        dataNameRE = fullfile(DRE(1).folder, DRE(1).name);
        ccep = load(dataNameRE);


        dataBase(i).sub_label = sub_label;
        dataBase(i).ses_label = ses_label;
        dataBase(i).tb_electrodes = tb_electrodes;
        dataBase(i).dwi = dwi;
        dataBase(i).ccep(j).run_label = run_label{j};
        dataBase(i).ccep(j).ccep = ccep;
        fprintf('...Subject %s %s has been run...\n',sub_label,run_label{j})
    end
end

disp('All subjects are loaded')

