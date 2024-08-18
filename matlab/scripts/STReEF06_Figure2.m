 %STReEF06_figure2
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 

% author: Susanne Jelsma
% date: October 2021

% % visualize the inter_modal similarity (figure 2 from the manuscript) 
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

%% Calculate the inter-modal similarity

% calculate the Jaccard Index
JI = zeros(size(dataBase,2,1)); JI_expected = zeros(size(dataBase,2,1));
for subj=1:size(dataBase,2)

[JI(subj), JI_expected(subj)] = jaccard(dataBase(subj).SC_matrix,dataBase(subj).EC_matrix);

end
clear subj 

% Calculate the ratio between structural and effective connections in the set of symmetric difference connections
ratio = zeros(size(dataBase,2,1)); 
for subj=1:size(dataBase,2)

ratio(subj) = ratio_SC_EC(dataBase(subj).SC_matrix,dataBase(subj).EC_matrix);

end

disp('Inter-modal similarity calculated')

%% visualize the inter_modal similarity (figure 2 from the manuscript) 

[R,J] = intermodal_similarity(JI, JI_expected, ratio);
saveas(R,'ratio','epsc') % save the figure for further processing with Adobe Illustrator
saveas(J,'Jaccard Index','epsc') % save the figure for further processing with Adobe Illustrator

clear ratio subj JI JI_expected R J