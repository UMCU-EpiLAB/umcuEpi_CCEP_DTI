 %STReEF06_figureS4
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 

% author: Susanne Jelsma
% date: October 2021

% % visualize the network topology for all patients volume of electrode contact areas (VEA) (figure S4 in manuscript)
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

%% load volume of electrode contact areas and BIDS electrodes information (electrodes.tsv)

dataBase = load_network_data(myDataPath,cfg);

disp('Data loaded')

%%  calculate the network topology measures

for subj = 1:size(dataBase,2)

dataBase(subj).topology = calculate_topology(dataBase(subj).SC_matrix,dataBase(subj).EC_matrix,dataBase(subj).tb_electrodes,dataBase(subj).elec_include);

end
disp('network topology calculated')

%% %  visualize the network topology for all patients volume of electrode contact areas (VEA) (figure S4 in manuscript)
i = 0;
for subj = 1%:5 % grid % plot in 3 parts to get the right dimensions

i = i + 1;
% compare with degree structural connectivity
VS = visual_topology_predictor(dataBase(subj).VEA, dataBase(subj).topology.degree_SC, dataBase(subj).tb_electrodes,dataBase(subj).elec_include, i);

end
saveas(VS,'correlation VEA structural grid','epsc') % save the figure for further processing with Adobe Illustrator

i = 0;
for subj = 6:11 %seeg part 1  % plot in 3 parts to get the right dimensions

i = i + 1;
% compare with degree structural connectivity
VS = visual_topology_predictor(dataBase(subj).VEA, dataBase(subj).topology.degree_SC, dataBase(subj).tb_electrodes,dataBase(subj).elec_include, i);

end
saveas(VS,'correlation VEA structural seeg 1','epsc') % save the figure for further processing with Adobe Illustrator

i = 0;
for subj = 12:13 %seeg part 2 % plot in 3 parts to get the right dimensions

i = i + 1;
% compare with degree structural connectivity
VS = visual_topology_predictor(dataBase(subj).VEA, dataBase(subj).topology.degree_SC, dataBase(subj).tb_electrodes,dataBase(subj).elec_include, i);

end
saveas(VS,'correlation VEA structural seeg 2','epsc') % save the figure for further processing with Adobe Illustrator

