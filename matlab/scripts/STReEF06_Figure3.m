 %STReEF06_figure3
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 
% detect CCEPs in either cECoG or sEEG data

% author: Susanne Jelsma
% date: October 2021

% % visualize the inter_modal similarity with connectivity matrices (figure 3 from the manuscript) 
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

%% detected and visual checked CCEP data and BIDS electrodes information (electrodes.tsv)

dataBase = load_network_data(myDataPath,cfg);

disp('Detected data loaded')

%% calculate the Jaccard Index
JI = zeros(size(dataBase,2,1)); JI_expected = zeros(size(dataBase,2,1));
for subj=1:size(dataBase,2)

[JI(subj), JI_expected(subj)] = jaccard(dataBase(subj).SC_matrix,dataBase(subj).EC_matrix);

end
clear subj 

disp('Jaccard Index calculated')

%% visualize the inter-modal similarity with connectivity matrices (figure 3 from the manuscript)

for subj = [5,4] %size(dataBase,2)

V = visual_networks(dataBase(subj).SC_matrix,dataBase(subj).EC_matrix,dataBase(subj).sub_label); % function to plot the connectivity matrices
V.WindowState = 'maximized';
saveas(V,sprintf('visual_networks%s',dataBase(subj).sub_label),'epsc') % save the figure for further processing with Adobe Illustrator

end

clear subj V