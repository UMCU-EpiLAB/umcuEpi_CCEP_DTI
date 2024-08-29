 %STReEF06_figureS5
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 

% author: Susanne Jelsma
% date: October 2021

% % visualize the network topology for all patients betweenness centrality (figure S5 manuscript)
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

%% structural and effective connectivity

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
    dataBase(nSubj).modality = EC.modality;
    dataBase(nSubj).EC_matrix = EC.EC_matrix;
    dataBase(nSubj).SC_matrix = SC.SC_matrix;
    dataBase(nSubj).x_select = EC.x_select;
    dataBase(nSubj).y_select = EC.y_select;
    dataBase(nSubj).z_select = EC.z_select;
    dataBase(nSubj).VEA = SC.VEA;

end

% housekeeping
clear EC fileName fileName2 nSubj SC ses_label sub_label

disp('Data loaded')

%%  calculate the network topology measures

dataBase = calculate_topology(dataBase);

disp('Network topology calculated')

%% visualize the network topology for all patients betweenness centrality (figure S5 manuscript)

count = 0;
ecog = find(contains({dataBase(:).modality},'ecog'));
for nSubj =  ecog %grid  % plot in 3 parts to get the right dimensions

count = count + 1;
% compare with degree structural connectivity
T = visual_topology_predictor(dataBase(nSubj).topology.BC_EC, dataBase(nSubj).topology.BC_SC,dataBase(nSubj).soz_select,count);

end
saveas(T,'correlation BC grid','epsc') % save the figure for further processing with Adobe Illustrator

close all 

count = 0;
seeg = find(contains({dataBase(:).modality},'seeg'));
for nSubj = seeg(1:6) %seeg part 1 %plot in 3 parts to get the right dimensions

count= count+ 1;
% compare with degree structural connectivity
T = visual_topology_predictor(dataBase(nSubj).topology.BC_EC, dataBase(nSubj).topology.BC_SC, dataBase(nSubj).soz_select,count);

end
saveas(T,'correlation BC seeg 1','epsc') % save the figure for further processing with Adobe Illustrator

close all

count = 0;
for nSubj = seeg(7:end) % seeg part 2 % plot in 3 parts to get the right dimensions

count= count+ 1;
% compare with degree structural connectivity
T = visual_topology_predictor(dataBase(nSubj).topology.BC_EC, dataBase(nSubj).topology.BC_SC, dataBase(nSubj).soz_select,count);

end
saveas(T,'correlation BC seeg 2','epsc') % save the figure for further processing with Adobe Illustrator


