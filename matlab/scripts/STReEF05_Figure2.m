%STReEF05_figure2
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al.

% author: Susanne Jelsma
% date: October 2021

% This script:
% - set correct paths
% - define subjects to include
% - load preprocessed structural and effective matrices (STReEF02_postprocessEC.m)
% - calculate intermodal similarity
% - visualize the inter_modal similarity (figure 2 from the manuscript)

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

%% load structural and effective connectivity

dataBase = struct();

for nSubj = 1:size(cfg.sub_label,2)

    sub_label = ['sub-' cfg.sub_label{nSubj}];
    ses_label = cfg.ses_label{nSubj};

    fileName = [sub_label,'_',ses_label,'_Effective_Connectivity.mat'];
    fileName2 = [sub_label,'_',ses_label,'_Structural_Connectivity.mat'];
    EC = load(fullfile(myDataPath.input_dev,sub_label,fileName));
    SC = load(fullfile(myDataPath.input_dev,sub_label,fileName2));

    dataBase(nSubj).sub_label = sub_label;
    dataBase(nSubj).ses_label = ses_label;

    dataBase(nSubj).modality = EC.modality;
    dataBase(nSubj).EC_matrix = EC.EC_matrix;
    dataBase(nSubj).SC_matrix = SC.SC_matrix;
    
end

% housekeeping
clear EC fileName fileName2 nSubj SC ses_label sub_label

disp('Data loaded')

%% Calculate the inter-modal similarity

% calculate the Jaccard Index
JI = zeros(size(dataBase,2,1)); JI_expected = zeros(size(dataBase,2,1));

for nSubj = 1:size(dataBase,2)

    [JI(nSubj), JI_expected(nSubj)] = jaccard(dataBase(nSubj).SC_matrix,dataBase(nSubj).EC_matrix);

end

% housekeeping
clear nSubj

% Calculate the ratio between structural and effective connections in the set of symmetric difference connections
ratio = zeros(size(dataBase,2,1));
for nSubj = 1:size(dataBase,2)

    ratio(nSubj) = ratio_SC_EC(dataBase(nSubj).SC_matrix,dataBase(nSubj).EC_matrix);

end

disp('Inter-modal similarity calculated')

%% visualize the inter_modal similarity (figure 2 from the manuscript)

[R,J] = intermodal_similarity(JI, JI_expected, ratio, {dataBase.modality});

% Create the folder if it doesn't exist already.
targetFolder = fullfile(myDataPath.output,'/Figures/');
if ~exist(targetFolder, 'dir')
    mkdir(targetFolder);
end

saveas(R,fullfile(targetFolder,'ratio'),'epsc') % save the figure for further processing with Adobe Illustrator
saveas(J,fullfile(targetFolder,'Jaccard Index'),'epsc') % save the figure for further processing with Adobe Illustrator

clear ratio nSubj JI JI_expected R J