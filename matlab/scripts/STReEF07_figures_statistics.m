% STReEF07_tables
% make figures with data gathered from statistics in R
% author: Susanne Jelsma 
% date:  feb 2024

%linear mixed model
data_llm1 = [2.259;-3.695;0.091;-1.010]; %11 subj Degree Structural Connectivity ~ 1 * ß0+ Degree Effective Connectivity * ßDEC+ Node Proximity * ßNP + Volume Electrode Area * ßVEA + Seizure Onset Zone nodes * ßSOZ + (1 * ß0 | Subject)Degree Structural Connectivity ~ 1 * ß0+ Degree Effective Connectivity * ßDEC+ Node Proximity * ßNP + Volume Electrode Area * ßVEA + Seizure Onset Zone nodes * ßSOZ + (1 * ß0 | Subject)
data_llm2 = [2.258;-3.7;NaN;-1.010];
data_llm3 = [2.190;-3.585;NaN;NaN];
data_all = [2.425;-4.009;NaN;NaN];
stat_tresh = 1.647131;
stat_tresh_min = -1.647131;

color = [204/250 37/250 41/250]; %red
clbar = [128/250 133/250 133/250]; %grey

f1 = figure(1);
ylabel('Test statistic T linear multilevel model')

subplot(3,1,1)
bar(data_llm1,'w')
set(gca, 'XTickLabel',[{'DEC'} {'NP'} {'VEA'} {'SOZ'}])
ylim([-4,4])
set(gca, 'FontSize',12) %12
yline(stat_tresh, Color=color,LineWidth=1.5)
yline(stat_tresh_min, Color=color,LineWidth=1.5)

subplot(3,1,2)
bar(data_llm2,'w')
set(gca, 'XTickLabel',[{'DEC'} {'NP'} {'VEA'} {'SOZ'}])
ylim([-4,4])
set(gca, 'FontSize',12) %12
yline(stat_tresh, Color=color,LineWidth=1.5)
yline(stat_tresh_min, Color=color,LineWidth=1.5)

subplot(3,1,3)
bar(data_all,'w')
set(gca, 'XTickLabel',[{'DEC'} {'NP'} {'VEA'} {'SOZ'}])
ylim([-4,4])
set(gca, 'FontSize',12) %12
yline(stat_tresh, Color=color,LineWidth=1.5)
yline(stat_tresh_min, Color=color,LineWidth=1.5)

saveas(f1,'lmm2','epsc') % save the figure for further processing with Adobe Illustrator
clear clbar color data_all data_llm1 data_llm2 data_llm3 f1 stat_tresh stat_tresh_min
%%
%jaccard index
pt = 1:13;
JI_exp= [0.14;0.14;0.23;0.10;0.19;0.08;0.05;0.15;0.05;0.08;0.20;0.06;0.07];
JI_obsv = [0.19;0.26;0.38;0.18;0.37;0.27;0.17;0.28;0.22;0.19;0.30;0.23;0.25];

f1= figure(1);
s = scatter(pt,JI_exp,300,'LineWidth', 1.8); 
color = [0.99 0.38 0.22];
s.MarkerEdgeColor = color;
s.Marker = '_'; % markers(subj);

hold on
s2 = scatter(pt,JI_obsv,300,'LineWidth', 1.8); 
color = [0.99 0.38 0.22]; % orange
s2.MarkerEdgeColor = color;
s2.Marker = '.'; % markers(subj);


xlim([0,13])
ylim([0,1])

set(gca, 'XTick',1:13)%numelxt
set(gca, 'YTick',[0,0.5,1])%numelxt

set(gca, 'FontSize',16) %12

ylabel('Jaccard Index')
xlabel('Patient')
legend(' Expected Jaccard Index', ' Observed Jaccard Index','FontSize',16)
saveas(f1,'Jaccard Index','epsc') % save the figure for further processing with Adobe Illustrator

clear pt JI_exp JI_obsv color sw s f1 s2