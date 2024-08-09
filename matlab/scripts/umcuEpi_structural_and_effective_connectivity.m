
% umcuEpi_structural_and_effective_connectivity
% inter-modal similarity and network topology comparison between structural (derived from DWI) and effective (derived from SPES) networks

% author: Susanne Jelsma & Dorien van Blooijs
% date: Aug 2024

% This matlab code is developed for the manuscript 'Structural and
% Effective brain connectivity in focal epilepsy' by Jelsma et al. 

%% SECTION 1: load data 

% set paths

clc
clear
cfg.folderinput = 'DATA_for_publication'; % from which folder would you like to loadthe main data?
myDataPath = setLocalDataPath(cfg);

% load data

load([myDataPath.dataPath 'dataMain.mat']);
%% SECTION 2: Inter-modal similarity

% calculate the Jaccard Index
JI = zeros(size(dataMain,2,1)); JI_expected = zeros(size(dataMain,2,1));
for subj=1:size(dataMain,2)

[JI(subj), JI_expected(subj)] = jaccard_SJ(dataMain(subj).SC_matrix,dataMain(subj).EC_matrix);

end
clear subj 

% Calculate the ratio between structural and effective connections in the set of symmetric difference connections
ratio = zeros(size(dataMain,2,1)); 
for subj=1:size(dataMain,2)

ratio(subj) = ratio_SC_EC(dataMain(subj).SC_matrix,dataMain(subj).EC_matrix);

end

% visualize the inter_modal similarity (figure 2 from the manuscript) 
[R,J] = intermodal_similarity(JI, JI_expected, ratio);
saveas(R,'ratio','epsc') % save the figure for further processing with Adobe Illustrator
saveas(J,'Jaccard Index','epsc') % save the figure for further processing with Adobe Illustrator

% visualize the inter-modal similarity with connectivity matrices (figure 3
% from the manuscript)
for subj = 1:size(dataMain,2)

V = visual_networks_SJ(dataMain(subj).SC_matrix,dataMain(subj).EC_matrix,dataMain(subj).index); % function to plot the connectivity matrices
V.WindowState = 'maximized';
saveas(V,sprintf('visual_networks%d',dataMain(subj).index),'epsc') % save the figure for further processing with Adobe Illustrator

end

clear ratio subj JI JI_expected R J V

%% SECTION 3: Network topology

% calculate the network topology measures and node proximity
network_topology = struct('topology', cell(1, size(dataMain,2)));
for subj = 1:size(dataMain,2)

network_topology(subj).topology = calculate_topology(dataMain(subj).SC_matrix,dataMain(subj).EC_matrix,dataMain(subj).tb_electrodes,dataMain(subj).elec_include);

end

%  visualize the network topology for all patients degree (figure 4 manuscript)
i = 0;
for subj = 12:13 % 1:5 grid % 6:11 seeg part 1 % 12:13 seeg part 2 % plot in 3 parts to get the right dimensions

i = i + 1;
T = visual_topology(network_topology(subj).topology.degree_EC, network_topology(subj).topology.degree_SC, dataMain(subj).tb_electrodes,dataMain(subj).elec_include, i);

end
saveas(T,'correlation degree grid','epsc') % save the figure for further processing with Adobe Illustrator

%  visualize the network topology for all patients node proximity (figures S1 and S2 manuscript)

i = 0;
for subj = 1:5 % 1:5 grid % 6:11 seeg part 1 % 12:13 seeg part 2 % plot in 3 parts to get the right dimensions

i = i + 1;
% compare with degree structural connectivity
NS = visual_topology_predictor(network_topology(subj).topology.node_proximity, network_topology(subj).topology.degree_SC, dataMain(subj).tb_electrodes,dataMain(subj).elec_include, i);

end
saveas(NS,'correlation np structural grid','epsc') % save the figure for further processing with Adobe Illustrator

i = 0;
for subj = 1:5 % 1:5 grid % 6:11 seeg part 1 % 12:13 seeg part 2 % plot in 3 parts to get the right dimensions

i = i + 1;
% compare with degree effective connectivity
NE = visual_topology_predictor(network_topology(subj).topology.node_proximity, network_topology(subj).topology.degree_EC, dataMain(subj).tb_electrodes,dataMain(subj).elec_include, i);

end
saveas(NE,'correlation np effective grid','epsc') % save the figure for further processing with Adobe Illustrator

% Histogram of the volume of electrode contact areas (figure S3 manuscript)
i = 0;
for subj = 1:5 % 1:5 grid % 6:11 seeg part 1 % 12:13 seeg part 2 % plot in 3 parts to get the right dimensions

i = i + 1;
% compare with degree structural connectivity
VEA = visual_VEA(dataMain(subj).VEA, dataMain(subj).tb_electrodes,dataMain(subj).elec_include, i);

end
saveas(VEA,'histogram VEA grid','epsc') % save the figure for further processing with Adobe Illustrator

%  visualize the network topology for all patients volume of electrode contact areas (VEA) (figure S4 in manuscript)
i = 0;
for subj = 1:5 % 1:5 grid % 6:11 seeg part 1 % 12:13 seeg part 2 % plot in 3 parts to get the right dimensions

i = i + 1;
% compare with degree structural connectivity
VS = visual_topology_predictor(dataMain(subj).VEA, network_topology(subj).topology.degree_SC, dataMain(subj).tb_electrodes,dataMain(subj).elec_include, i);

end
saveas(VS,'correlation VEA structural seeg2','epsc') % save the figure for further processing with Adobe Illustrator

%  visualize the network topology for all patients betweenness centrality (figure S5 manuscript)

i = 0;
for subj = 12:13 % 1:5 grid % 6:11 seeg part 1 % 12:13 seeg part 2 % plot in 3 parts to get the right dimensions

i = i + 1;
T = visual_topology(network_topology(subj).topology.BC_EC, network_topology(subj).topology.BC_SC, dataMain(subj).tb_electrodes,dataMain(subj).elec_include, i);

end
saveas(T,'correlation BC grid','epsc') % save the figure for further processing with Adobe Illustrator
%% SECTION 4 : linear multilevel model
% prepare data for multilevel model calculations in R

% calculate the nr of channels over all patients
nr_channels =  NaN(size(dataMain,2),1);
for subj = 1:size(dataMain,2)
nr_channels(subj) = size(dataMain(subj).EC_matrix,1);
end
size_long = sum(nr_channels);

% make a matrix with a row for each channel
data_long = NaN(size_long,6); % matrix with for every channel per patient the predictors
i=1;
for subj = 1:size(dataMain,2)
soz_all = strcmpi(dataMain(subj).tb_electrodes.soz,'yes'); % electrodes in soz
soz = soz_all(dataMain(subj).elec_include); % included electrodes in soz

sz = size(dataMain(subj).EC_matrix,1); % nr pf channels

data_long(i:i+sz-1,1)= subj; % patient index
data_long(i:i+sz-1,2)= network_topology(subj).topology.degree_SC; % degree structural networks
data_long(i:i+sz-1,3)= network_topology(subj).topology.degree_EC; % degree effective networks
data_long(i:i+sz-1,4)= network_topology(subj).topology.node_proximity; % node proximity
data_long(i:i+sz-1,5)= soz;  %if the channel is in the SOZ or not
data_long(i:i+sz-1,6)= dataMain(subj).VEA; % volume of electrode contact areas

i=i+sz;
end


names_data = {'subj' , 'SCD' , 'ECD' , 'NP' , 'SOZ' , 'VEA'};
data_l=array2table(data_long,'VariableNames',names_data);

dataPath = myDataPath.DWIMATLABpath;  % change!
writetable(data_l,[dataPath 'data_long_LMM.csv']) % save the matrix with for every channel per patient the specifications

clear data_long data_l data_Path sz soz soz_all subj i names_data nr_channels size_long

% R warning
warning('run now R code in STReEF05_compare_networks_R and run section 2)')

% linear mixed model (verbeter nog!)
data_lmm = readtable("LLMdata_R.csv");
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


