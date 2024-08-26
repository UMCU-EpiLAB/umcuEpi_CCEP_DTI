 %STReEF06_figureS3
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 

% author: Susanne Jelsma
% date: October 2021

% % make histogram of the volume of electrode contact areas (figure S3 manuscript)
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

%% Histogram of the volume of electrode contact areas (figure S3 manuscript)

i = 0;
for subj = 1:5 % grid % plot in 3 parts to get the right dimensions

i = i + 1;
% compare with degree structural connectivity
VEA = visual_VEA(dataBase(subj).VEA, dataBase(subj).tb_electrodes,dataBase(subj).elec_include, i);

end
saveas(VEA,'histogram VEA grid','epsc') % save the figure for further processing with Adobe Illustrator

i = 0;
for subj = 6:11 %seeg part 1 % plot in 3 parts to get the right dimensions

i = i + 1;
% compare with degree structural connectivity
VEA = visual_VEA(dataBase(subj).VEA, dataBase(subj).tb_electrodes,dataBase(subj).elec_include, i);

end
saveas(VEA,'histogram VEA seeg 1','epsc') % save the figure for further processing with Adobe Illustrator

i = 0;
for subj = 12:13 %seeg part 2 % plot in 3 parts to get the right dimensions

i = i + 1;
% compare with degree structural connectivity
VEA = visual_VEA(dataBase(subj).VEA, dataBase(subj).tb_electrodes,dataBase(subj).elec_include, i);

end
saveas(VEA,'histogram VEA seeg 2','epsc') % save the figure for further processing with Adobe Illustrator
