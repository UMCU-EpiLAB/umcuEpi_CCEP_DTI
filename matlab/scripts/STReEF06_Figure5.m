 %STReEF06_figure5
% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 

% author: Susanne Jelsma
% date: October 2021

% % visualize results of linear mixed model (figure 5 from the manuscript) 

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

%% load linear mixed model data

data_lmm = readtable(fullfile(myDataPath.input_dev,'output_LMM_model.csv'));

disp('llm data loaded')

%% visualize results of linear mixed model
data_lmm1 = data_lmm.lmm1_t_value(2:5); % extract the data of the predictors and not the intercept (the first row)
data_lmm2 = data_lmm.lmm2_t_value(2:5); 
data_lmm3 = data_lmm.lmm3_t_value(2:5);
data_all = data_lmm.lmm_final_t_value(2:5);

stat_tresh = str2double(data_lmm.statistical_tresh(1)); 

color = [204/250 37/250 41/250]; %red

f1 = figure(1);

subplot(3,1,1)
bar(data_lmm1,'w')
set(gca, 'XTickLabel',[{'ECD'} {'NP'} {'SOZ'} {'VEA'}])
ylim([-4,4])
set(gca, 'FontSize',12) 
yline(stat_tresh, Color=color,LineWidth=1.5)
yline(-stat_tresh, Color=color,LineWidth=1.5)
ylabel('T-value')

subplot(3,1,2)
bar(data_lmm2,'w')
set(gca, 'XTickLabel',[{'ECD'} {'NP'} {'SOZ'} {'VEA'}])
ylim([-4,4])
set(gca, 'FontSize',12) 
yline(stat_tresh, Color=color,LineWidth=1.5)
yline(-stat_tresh, Color=color,LineWidth=1.5)
ylabel('T-value')

subplot(3,1,3)
bar(data_lmm3,'w')
set(gca, 'XTickLabel',[{'ECD'} {'NP'} {'SOZ'} {'VEA'}])
ylim([-4,4])
set(gca, 'FontSize',12) 
yline(stat_tresh, Color=color,LineWidth=1.5)
yline(-stat_tresh, Color=color,LineWidth=1.5)
ylabel('T-value')

f2 = figure(2);

subplot(3,1,1)
bar(data_all,'w')
set(gca, 'XTickLabel',[{'ECD'} {'NP'} {'SOZ'} {'VEA'}])
ylim([-4,4])
set(gca, 'FontSize',12) 
yline(stat_tresh, Color=color,LineWidth=1.5)
yline(-stat_tresh, Color=color,LineWidth=1.5)
ylabel('T-value')

saveas(f1,'lmm_figure5a','epsc') % save the figure for further processing with Adobe Illustrator
saveas(f2,'lmm_figure5b','epsc') % save the figure for further processing with Adobe Illustrator
clear color data_all data_lmm1 data_lmm2 data_lmm3 f1 f2 stat_tresh 

