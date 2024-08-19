 %STReEF03_postprocessEC
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 

% author: Susanne Jelsma
% date: October 2021

% load detected and visual checked CCEP data, merge the BIDS runs, include electrodes, make a symmetric bi-directional effective connectivity network
%% set paths
% set umcuEpi_CCEP_DTI/matlab in your directory and run this section

clc
clear

rootPath = matlab.desktop.editor.getActiveFilename;
RepoPath = fileparts(rootPath);
matlabFolder = strfind(RepoPath,'matlab');
addpath(genpath(RepoPath(1:matlabFolder+6)));

cfg.folderinput = 'shareData_STReEF'; % from which folder would you like to load ECoGs?
myDataPath = setLocalDataPath(cfg);

% housekeeping 
clear rootPath RepoPath matlabFolder

%% patient characteristics
sub_label = input('Patient number (STREEFXX) (select multiple patients by separating each with a comma): ','s');

cfg.sub_label = strsplit(sub_label,{', ',','});

cfg = selectPatients(cfg, myDataPath);

%% detected and visual checked CCEP data and BIDS electrodes information (electrodes.tsv and channels.tsv)

dataBase = load_derivatives_data(myDataPath,cfg);

disp('Detected data loaded')

%% merge runs

for subj= 1:size(dataBase,2)
    if size(cfg.run_label{subj},2) > 1
                dataBase(subj).ccep_m = dataBase(subj).ccep; % copy everything to ccep_merge to keep the data
                data = dataBase(subj).ccep(1); % start with the basic info from the first run
                dataBase(subj).ccep = data;
              % append the vectors cc_stimsets (containing channel number of stimulated electrode pairs) and cc_stimchans (containing channel names of stimulated electrode pairs) from each run into one large vector
            for run = 2:size(cfg.run_label{subj},2)
                n1_peak_sample = dataBase(subj).ccep_m(run).ccep.n1_peak_sample;
                n1_peak_amplitude = dataBase(subj).ccep_m(run).ccep.n1_peak_amplitude;
                scored = dataBase(subj).ccep_m(run).ccep.checked;
                stimsets = dataBase(subj).ccep_m(run).ccep.cc_stimsets;
                stimchannels = dataBase(subj).ccep_m(run).ccep.cc_stimchans;


                dataBase(subj).ccep.ccep.cc_stimsets = [dataBase(subj).ccep.ccep.cc_stimsets;stimsets];
                dataBase(subj).ccep.ccep.cc_stimchans = [dataBase(subj).ccep.ccep.cc_stimchans;stimchannels];                                         
                dataBase(subj).ccep.ccep.checked = [dataBase(subj).ccep.ccep.checked scored];
                dataBase(subj).ccep.ccep.n1_peak_sample = [dataBase(subj).ccep.ccep.n1_peak_sample n1_peak_sample];
                dataBase(subj).ccep.ccep.n1_peak_amplitude = [dataBase(subj).ccep.ccep.n1_peak_amplitude n1_peak_amplitude];
            
                % check for equal stim par
                run1 = dataBase(subj).ccep_m(1).ccep; run2 = dataBase(subj).ccep_m(run).ccep;             
                fields = {'n1_peak_sample','n1_peak_amplitude','checked','cc_stimsets','cc_stimchans','checkUntilStimp','dataName'};
                run1 = rmfield(run1,fields);run2 = rmfield(run2,fields);
                if ~isequal(dataBase(subj).ccep_m(1).ccep.ch,dataBase(subj).ccep_m(run).ccep.ch)
                        warning('channels %s are not equal to the first run',dataBase(subj).ccep(run).run_label)
                end            
                if ~isequal(run1,run2)
                        warning('parameters %s are not equal to the first run',dataBase(subj).ccep_m(run).run_label)
                end
            end
        fprintf('...Runs subject %s has been merged ... \n',dataBase(subj).sub_label)
    end
end

disp('runs merged')
clear i data run subj stimsets stimchannels stimchannels run1 run2 fields scored n1_peak_amplitude n1_peak_sample sub_label 

%% select and include the electrode channels

for  subj=1:size(dataBase,2)
% if sEEG, use only the gray matter channels, hippocampus, amygdala, lesion, and gliosis
if any(contains(fieldnames(dataBase(subj).tb_electrodes),'graymatter')) % select sEEG patients
    idx_screw = strcmpi(dataBase(subj).tb_electrodes.screw,'yes'); % remove screw electrodes
    idx_csf = strcmpi(dataBase(subj).tb_electrodes.csf,'yes'); % remove csf electrodes
    idx_whitematter = strcmpi(dataBase(subj).tb_electrodes.whitematter,'yes') & strcmpi(dataBase(subj).tb_electrodes.graymatter,'no'); % remove white matter channels (not the bordeline gray/white matter)
    idx_bad = strcmpi(dataBase(subj).ccep(1).tb_channels.status,'bad'); % remove bad channels from analysis because you cannot make extract CCEPs from it, so cannot compare it to the structural networks
    idx_all = sum([idx_screw, idx_csf, idx_whitematter, idx_bad],2);

    elec_include = cell(size(dataBase(subj).tb_electrodes,1),1);
    [elec_include{idx_all>0}] = deal('no'); % no in elec_include vector if it is a screw, csf, whitematter, or bad channel
    [elec_include{idx_all==0}] = deal('yes'); % yes in elec_include vector if it is something else (options: graymatter, hippocampus, amygdala, lesion, or gliosis)
elseif any(contains(dataBase(subj).tb_electrodes.group,'grid')) % select grid patients
    idx_silicon = strcmpi(dataBase(subj).tb_electrodes.silicon,'yes');% remove electrodes who are laying on other electrodes (silicon)
    idx_all = sum(idx_silicon,2);


    elec_include = cell(size(dataBase(subj).tb_electrodes,1),1);
    [elec_include{idx_all>0}] = deal('no');
    [elec_include{idx_all==0}] = deal('yes');
else
    warning('%s is not a sEEG or grid',dataBase(subj).sub_label)
    elec_include = cell(size(dataBase(subj).tb_electrodes,1),1);
    [elec_include{1:end}] = deal('yes');
end

% save in the dataBase for further use
elec_indx = strcmpi(elec_include,'yes');
dataBase(subj).ccep.elec_include = elec_indx;% save a logical vector including all channels and their inclusion yes (1) or no (0) in elec_include for easy computational matters

end

clear stimchannels chan_include elec_include elec_indx idx_all idx_bad bad idx_csf idx_screw idx_silicon idx_whitematter j i subjs run

%% make a symmetric bi-directional effective connectivity network
dataBase = construct_network(dataBase);

disp('networks constructed')

%% save effective network

for subj = 1:size(dataBase,2)
    targetFolder = fullfile(myDataPath.output, dataBase(subj).sub_label);

    % Create the folder if it doesn't exist already.
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end
    
    
    fileName = [dataBase(subj).sub_label,'_',dataBase(subj).ses_label,'_Effective_Connectivity_new.mat']; 

    % insert variables
    EC.elec_include = dataBase(subj).ccep.elec_include;
    EC.EC_matrix = dataBase(subj).ccep.EC_matrix;

    save(fullfile(targetFolder,fileName), '-struct','EC')
    fprintf('Saved network-struct in %s \n',fullfile(targetFolder,fileName))
end
clear fileName targetFolder EC_matrix