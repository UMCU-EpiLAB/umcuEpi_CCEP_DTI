%STReEF03_paramModel

% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al.

% author: Susanne Jelsma
% date: October 2021

% This script:
% - set correct paths
% - define subjects to include
% - load preprocessed CCEP data (STReEF01_preprocessCCEP.m)
% - determine network topology
% - prepare data for multilevel model calculations in R
% - save as csv-file

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

    dataBase(nSubj).ch_select = EC.ch_select;
    dataBase(nSubj).soz_select = EC.soz_select;
    dataBase(nSubj).EC_matrix = EC.EC_matrix;
    dataBase(nSubj).SC_matrix = SC.SC_matrix;
    dataBase(nSubj).x_select = EC.x_select;
    dataBase(nSubj).y_select = EC.y_select;
    dataBase(nSubj).z_select = EC.z_select;
    dataBase(nSubj).VEA = SC.VEA;

end

% housekeeping
clear EC fileName fileName2 nSubj SC ses_label sub_label

disp('Detected data loaded')

%% Network topology

% calculate the network topology measures and node proximity
dataBase = calculate_topology(dataBase);

disp('Network topology calculated')

%% prepare data for multilevel model calculations in R

% calculate the nr of channels for all patients combined
nr_channels =  NaN(size(dataBase,2),1);
for nSubj = 1:size(dataBase,2)
    nr_channels(nSubj) = size(dataBase(nSubj).EC_matrix,1);
end
size_long = sum(nr_channels);

% make a matrix with a row for each channel
data_long = NaN(size_long,6); % matrix with for every channel per patient the predictors

count = 1;
for nSubj = 1:size(dataBase,2)
    sizeSubj = size(dataBase(nSubj).EC_matrix,1); % nr pf channels

    data_long(count:count+sizeSubj-1,1)= nSubj; % patient index
    data_long(count:count+sizeSubj-1,2)= dataBase(nSubj).topology.degree_SC; % degree structural networks
    data_long(count:count+sizeSubj-1,3)= dataBase(nSubj).topology.degree_EC; % degree effective networks
    data_long(count:count+sizeSubj-1,4)= dataBase(nSubj).topology.node_proximity; % node proximity
    data_long(count:count+sizeSubj-1,5)= dataBase(nSubj).soz_select;  %if the channel is in the SOZ or not
    data_long(count:count+sizeSubj-1,6)= dataBase(nSubj).VEA; % volume of electrode contact areas

    count = count+sizeSubj;
end

names_data = {'subj' , 'SCD' , 'ECD' , 'NP' , 'SOZ' , 'VEA'};
tb_data_long = array2table(data_long,'VariableNames',names_data);

writetable(tb_data_long,[myDataPath.output 'input_LMM_model_new.csv']) % save the matrix with for every channel per patient the specifications

% housekeeping
clear count data_long tb_data_long sizeSubj nSubj count names_data nr_channels size_long
