 %STReEF04_paramModel
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 

% author: Susanne Jelsma
% date: October 2021

%  prepare data for multilevel model calculations in R
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

dataBase = load_network_data(myDataPath,cfg);

disp('Detected data loaded')

%% Network topology

% calculate the network topology measures and node proximity
for subj = 1:size(dataBase,2)

dataBase(subj).topology = calculate_topology(dataBase(subj).SC_matrix,dataBase(subj).EC_matrix,dataBase(subj).tb_electrodes,dataBase(subj).elec_include);

end
disp('network topology calculated')


%% prepare data for multilevel model calculations in R

% calculate the nr of channels over all patients
nr_channels =  NaN(size(dataBase,2),1);
for subj = 1:size(dataBase,2)
nr_channels(subj) = size(dataBase(subj).EC_matrix,1);
end
size_long = sum(nr_channels);

% make a matrix with a row for each channel
data_long = NaN(size_long,6); % matrix with for every channel per patient the predictors
i=1;
for subj = 1:size(dataBase,2)
soz_all = strcmpi(dataBase(subj).tb_electrodes.soz,'yes'); % electrodes in soz
soz = soz_all(dataBase(subj).elec_include); % included electrodes in soz

sz = size(dataBase(subj).EC_matrix,1); % nr pf channels

data_long(i:i+sz-1,1)= subj; % patient index
data_long(i:i+sz-1,2)= dataBase(subj).topology.degree_SC; % degree structural networks
data_long(i:i+sz-1,3)= dataBase(subj).topology.degree_EC; % degree effective networks
data_long(i:i+sz-1,4)= dataBase(subj).topology.node_proximity; % node proximity
data_long(i:i+sz-1,5)= soz;  %if the channel is in the SOZ or not
data_long(i:i+sz-1,6)= dataBase(subj).VEA; % volume of electrode contact areas

i=i+sz;
end


names_data = {'subj' , 'SCD' , 'ECD' , 'NP' , 'SOZ' , 'VEA'};
data_l=array2table(data_long,'VariableNames',names_data);

dataPath = myDataPath.output;
writetable(data_l,[dataPath 'input_LMM_model_new.csv']) % save the matrix with for every channel per patient the specifications

clear data_long data_l data_Path sz soz soz_all subj i names_data nr_channels size_long
