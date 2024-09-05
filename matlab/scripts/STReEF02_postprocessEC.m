 %STReEF02_postprocessEC
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 

% author: Susanne Jelsma
% date: October 2021

% This script is used to:
% - set correct paths
% - define subjects to include
% - load detected and visual checked CCEP data
% - merge multiple SPES-runs
% - make a symmetric bi-directional effective connectivity network
% - save matrix in derivatives

%% set paths
% set umcuEpi_CCEP_DTI/matlab in your directory and run this section

clc
clear

rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

myDataPath = STReEF_setLocalDataPath(1);

% housekeeping 
clear rootPath RepoPath matlabFolder

%% patient characteristics
sub_label = ['STREEF01, STREEF02, STREEF03, STREEF04, STREEF05, STREEF06,' ...
    'STREEF07, STREEF08, STREEF09, STREEF10, STREEF11, STREEF12, STREEF13'];

cfg.sub_label = strsplit(sub_label,{', ',','});

cfg = selectPatients(cfg, myDataPath);

% housekeeping
clear sub_label

%% detected and visual checked CCEP data and BIDS electrodes information (electrodes.tsv and channels.tsv)

dataBase = load_derivatives_data(myDataPath,cfg);

disp('Data with detected and visually checked CCEP-N1s loaded')

%% merge runs

for nSubj = 1:size(dataBase,2)

    dataBase(nSubj).ccep.ch = dataBase(nSubj).metadata(1).ccep.ch;
    dataBase(nSubj).ccep.cc_stimsets = dataBase(nSubj).metadata(1).ccep.cc_stimsets;
    dataBase(nSubj).ccep.cc_stimchans = dataBase(nSubj).metadata(1).ccep.cc_stimchans;

    dataBase(nSubj).ccep.epoch_length = dataBase(nSubj).metadata(1).ccep.epoch_length;
    dataBase(nSubj).ccep.epoch_prestim = dataBase(nSubj).metadata(1).ccep.epoch_prestim;
    dataBase(nSubj).ccep.dir = dataBase(nSubj).metadata(1).ccep.dir;
    dataBase(nSubj).ccep.amp = dataBase(nSubj).metadata(1).ccep.amp;
    dataBase(nSubj).ccep.reref = dataBase(nSubj).metadata(1).ccep.reref;

    dataBase(nSubj).ccep.amplitude_thresh = dataBase(nSubj).metadata(1).ccep.amplitude_thresh;
    dataBase(nSubj).ccep.n1_peak_range = dataBase(nSubj).metadata(1).ccep.n1_peak_range;

    dataBase(nSubj).ccep.n1_peak_sample = dataBase(nSubj).metadata(1).ccep.n1_peak_sample;
    dataBase(nSubj).ccep.n1_peak_amplitude = dataBase(nSubj).metadata(1).ccep.n1_peak_amplitude;
    dataBase(nSubj).ccep.checked = dataBase(nSubj).metadata(1).ccep.checked;

    % append the vectors cc_stimsets (containing channel number of stimulated electrode pairs) and cc_stimchans (containing channel names of stimulated electrode pairs) from each run into one large vector
    for nRun = 2:size(cfg.run_label{nSubj},2)

        % check for equal settings in both runs:
        run1 = dataBase(nSubj).metadata(1).ccep; runN = dataBase(nSubj).metadata(nRun).ccep;
        fields = {'n1_peak_sample','n1_peak_amplitude','checked','cc_stimsets','cc_stimchans','checkUntilStimp','dataName'};
        run1 = rmfield(run1,fields); runN = rmfield(runN,fields);
        if ~isequal(run1,runN)
            error('Parameters %s are not equal to the first run: %s',...
                dataBase(nSubj).metadata(nRun).run_label,...
                dataBase(nSubj).metadata(1).run_label)
        else

            dataBase(nSubj).ccep.cc_stimsets = [dataBase(nSubj).ccep.cc_stimsets;...
                dataBase(nSubj).metadata(nRun).ccep.cc_stimsets;];
            dataBase(nSubj).ccep.cc_stimchans = [dataBase(nSubj).ccep.cc_stimchans;...
                dataBase(nSubj).metadata(nRun).ccep.cc_stimchans];
            dataBase(nSubj).ccep.checked = [dataBase(nSubj).ccep.checked ...
                dataBase(nSubj).metadata(nRun).ccep.checked];
            dataBase(nSubj).ccep.n1_peak_sample = [dataBase(nSubj).ccep.n1_peak_sample ...
                dataBase(nSubj).metadata(nRun).ccep.n1_peak_sample];
            dataBase(nSubj).ccep.n1_peak_amplitude = [dataBase(nSubj).ccep.n1_peak_amplitude ...
                dataBase(nSubj).metadata(nRun).ccep.n1_peak_amplitude];

        end
    end
    fprintf('... Runs %s has been merged ... \n',dataBase(nSubj).sub_label)
end

disp('All runs are merged')

% housekeeping
clear i data nRun nSubj stimsets stimchannels stimchannels run1 runN fields scored n1_peak_amplitude n1_peak_sample sub_label 


%% make a symmetric bi-directional effective connectivity network
dataBase = construct_connectivity(dataBase);

disp('Symmetrical, bidirectional connectivity matrices constructed')

%% select and include electrode channels that are located in/on graymatter, and not noisy upon visual evaluation

for  nSubj = 1:size(dataBase,2)

    % if sEEG, use only the gray matter channels, hippocampus, amygdala, lesion, and gliosis
    if any(contains(fieldnames(dataBase(nSubj).tb_electrodes),'graymatter')) % select sEEG patients
        idx_screw = strcmpi(dataBase(nSubj).tb_electrodes.screw,'yes'); % remove screw electrodes
        idx_csf = strcmpi(dataBase(nSubj).tb_electrodes.csf,'yes'); % remove csf electrodes
        idx_whitematter = strcmpi(dataBase(nSubj).tb_electrodes.whitematter,'yes') & strcmpi(dataBase(nSubj).tb_electrodes.graymatter,'no'); % remove white matter channels (not the bordeline gray/white matter)

        idx_bad = false(size(idx_csf,1),size(dataBase(nSubj).metadata,2));
        for nRun = 1:size(dataBase(nSubj).metadata,2)
            idx_bad(:,nRun) = strcmpi(dataBase(nSubj).metadata(nRun).tb_channels.status,'bad'); % remove bad channels from analysis because you cannot make extract CCEPs from it, so cannot compare it to the structural networks
        end

        idx_excl = sum([idx_screw, idx_csf, idx_whitematter, idx_bad],2);

        dataBase(nSubj).modality = 'seeg';

    elseif any(contains(dataBase(nSubj).tb_electrodes.group,'grid')) % select grid patients
        idx_silicon = strcmpi(dataBase(nSubj).tb_electrodes.silicon,'yes');% remove electrodes who are laying on other electrodes (silicon)

        idx_bad = false(size(idx_silicon,1),size(dataBase(nSubj).metadata,2));
        for nRun = 1:size(dataBase(nSubj).metadata,2)
            idx_bad(:,nRun) = strcmpi(dataBase(nSubj).metadata(nRun).tb_channels.status,'bad'); % remove bad channels from analysis because you cannot make extract CCEPs from it, so cannot compare it to the structural networks
        end

        idx_excl = sum([idx_silicon, idx_bad],2);
        
        dataBase(nSubj).modality = 'ecog';
    end

    elec_incl = false(size(idx_excl));
    elec_incl(idx_excl==0) = true;

    ch_select = dataBase(nSubj).ccep.ch(elec_incl);
    
    if isequal(dataBase(nSubj).tb_electrodes.name,dataBase(nSubj).ccep.ch)
        if iscell(dataBase(nSubj).tb_electrodes.x)
            dataBase(nSubj).ccep.x_select = str2double(dataBase(nSubj).tb_electrodes.x(elec_incl));
            dataBase(nSubj).ccep.y_select = str2double(dataBase(nSubj).tb_electrodes.y(elec_incl));
            dataBase(nSubj).ccep.z_select = str2double(dataBase(nSubj).tb_electrodes.z(elec_incl));
        else

            dataBase(nSubj).ccep.x_select = dataBase(nSubj).tb_electrodes.x(elec_incl);
            dataBase(nSubj).ccep.y_select = dataBase(nSubj).tb_electrodes.y(elec_incl);
            dataBase(nSubj).ccep.z_select = dataBase(nSubj).tb_electrodes.z(elec_incl);
        end
    else
        error('Channels names is unequal to names in tb_electrodes: %s',dataBase(nSubj).sub_label)
    end

    soz = strcmpi(dataBase(nSubj).tb_electrodes.soz,'yes');
    soz_select = soz(elec_incl);

    % remove all electrodes in the connectivity matrix if not in/on gray matter, of bad upon visual evaluation
    EC_matrix = dataBase(nSubj).ccep.EC_matrix;
    EC_matrix_select = EC_matrix(elec_incl,elec_incl);

    % save in the dataBase for further use
    dataBase(nSubj).ccep.ch_select = ch_select;
    dataBase(nSubj).ccep.soz_select = soz_select;
    dataBase(nSubj).ccep.EC_matrix_select = EC_matrix_select;
end

% housekeeping
clear EC_matrix EC_matrix_select elec_incl idx_excl idx_bad idx_csf idx_screw idx_silicon idx_whitematter nRun nSubj

%% save effective network

for nSubj = 1:size(dataBase,2)
    targetFolder = fullfile(myDataPath.output, dataBase(nSubj).sub_label);

    % Create the folder if it doesn't exist already.
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end
    
    fileName = [dataBase(nSubj).sub_label,'_',dataBase(nSubj).ses_label,'_Effective_Connectivity.mat'];

    % insert variables
    % EC.elec_incl = dataBase(nSubj).ccep.elec_incl;
    EC.EC_matrix = dataBase(nSubj).ccep.EC_matrix_select;
    EC.modality = dataBase(nSubj).modality;
    EC.ch_select = dataBase(nSubj).ccep.ch_select;
    EC.soz_select = dataBase(nSubj).ccep.soz_select;
    EC.x_select = dataBase(nSubj).ccep.x_select;
    EC.y_select = dataBase(nSubj).ccep.y_select;
    EC.z_select = dataBase(nSubj).ccep.z_select;

    save(fullfile(targetFolder,fileName), '-struct','EC')
    fprintf('Saved network-struct in %s \n',fullfile(targetFolder,fileName))
end

% housekeeping
clear fileName targetFolder EC nSubj