% STReEF05_compare_networks
% inter-modal similarity and network topography comparison between structural (derived from DWI) and effective (derived from SPES) networks

% author: Susanne Jelsma & Dorien van Blooijs
% date: May 2022

% Be aware! This file has a twin written in R code and to execute the code correctly, sections of the file: ' STReEF05_compare_networks_R' must be runned when indicated.

% load processed electrode contact area data, scored SPES/CCEP data, and SOZ information, construct an effective network matrix, make the effective network matrix symmetric, construct an structural network matrix,
% define the SOZ channels, calculate the inter-modal similarity, prepare for the jaccard index, calculate the network topography, calculate the node proximity, prepare data for multilevel model

%% SECTION 1: prepare structural and effective network data 

%% set paths
% set umcuEpi_CCEP_DTI/matlab in your directory and run this section

clc
clear
cfg.folderinput = 'chronic_ECoG'; % from which folder would you like to load ECoGs?
myDataPath = setLocalDataPath(cfg);

%% patient characteristics
%search for the SPES patients who have a DWI scan

files = dir(myDataPath.DWIpath);
files(1:2) = [];
sub_label = cell(1,size(files,1));
for subj=1:size(files,1)
    sub_label(1,subj)= cellstr(erase(files(subj).name,'sub-'));
end

i = 1;
for subj= 1:size(sub_label,2)
    x = input(sprintf('Load subject %s? (y/n): ',char(sub_label(subj))),'s');       
    if strcmp(x,'y') 
        cfg.sub_label(1,i) = sub_label(subj) ;
        i=i+1;
    end
end

clear sub_label files i x subj

cfg = selectPatients(cfg, myDataPath);

%% load processed electrode contact area data (see STReEF02_coreg_roidef_mrtrix) and visual scored CCEP data (see umcuEpi_CCEP_DTI (main)/matlab/scripts/ccepDTI03_detectCCEP.m for scoring CCEP data)

dataBase = load_derivatives_data(myDataPath,cfg);

disp('visual checked data loaded')

%% extract SOZ information
for subj = 1:size(cfg.sub_label,2)
dataBase(subj).soz.soz = strcmpi(dataBase(subj).tb_electrodes.soz,'yes');
end
disp('SOZ information loaded')

%% merge runs
% Be aware! Code with potential of causing errors due to copying of information, so check for new types of patients or new code added!
for subj= 1:size(dataBase,2)
    if size(dataBase(subj).ccep,2) > 1
       dataBase(subj).ccep_runs = dataBase(subj).ccep; % copy everything to ccep_runs
       data = dataBase(subj).ccep(1).ccep; % start with the basic info from the first run
       dataBase(subj).ccep = data;
       dataBase(subj).ccep.run_label = dataBase(subj).ccep_runs(1).run_label;
            for run = 2:size(cfg.run_label{subj},2)
                n1_peak_sample = dataBase(subj).ccep_runs(run).ccep.n1_peak_sample;
                n1_peak_amplitude = dataBase(subj).ccep_runs(run).ccep.n1_peak_amplitude;
                scored = dataBase(subj).ccep_runs(run).ccep.checked;
                stimsets = dataBase(subj).ccep_runs(run).ccep.cc_stimsets;
                stimchannels = dataBase(subj).ccep_runs(run).ccep.cc_stimchans;
                run_label = dataBase(subj).ccep_runs(run).run_label;                  
    
                dataBase(subj).ccep.run_label = [dataBase(subj).ccep.run_label run_label];
                dataBase(subj).ccep.checked = [dataBase(subj).ccep.checked scored];
                dataBase(subj).ccep.n1_peak_sample = [dataBase(subj).ccep.n1_peak_sample n1_peak_sample];
                dataBase(subj).ccep.n1_peak_amplitude = [dataBase(subj).ccep.n1_peak_amplitude n1_peak_amplitude];
                dataBase(subj).ccep.cc_stimsets = [dataBase(subj).ccep.cc_stimsets;stimsets];
                dataBase(subj).ccep.cc_stimchans = [dataBase(subj).ccep.cc_stimchans;stimchannels];
            end
        fprintf('...Runs subject %s has been merged ... \n',dataBase(subj).sub_label)
    else
       data = dataBase(subj).ccep.ccep; % if subject contains only one run
       run_label = dataBase(subj).ccep.run_label;
       dataBase(subj).ccep = data;
       dataBase(subj).ccep.run_label = run_label;
    end
end

disp('runs merged')
clear scored run subj stimsets stimchannels  n1_peak_sample n1_peak_amplitude data run_label

 %% SECTION 2: construct the effective and structural networks
dataBase = construct_network(dataBase,myDataPath);

disp('networks constructed')

%% SECTION 3: calculate the inter-modal similarity, prepare for the Jaccard Index calculation in R, calculate the ratio of structural and effective connections that did not overlap
JI = zeros(size(dataBase,2,1)); JI_expected = zeros(size(dataBase,2,1));JI_parker = zeros(size(dataBase,2,1));JI_alternative = zeros(size(dataBase,2,1));
for subj=1:size(dataBase,2)
dataPath = myDataPath.DWIMATLABpath; 
sub_label = erase(dataBase(subj).sub_label,'sub-RESP');
dir = [dataPath dataBase(subj).sub_label];

[JI(subj) JI_expected(subj) JI_parker(subj) JI_alternative(subj)] = jaccard(dataBase(subj).dwi.SC_matrix,dataBase(subj).ccep.EC_matrix,dir);

end
clear dataPath sub_label dir subj 

% ratio of structural and effective connections that did not overlap
% The ratio between structural and effective connections in the set of symmetric difference connections was calculated with: 
% Ratio_SC∆EC=log_10⁡((SC-(SC⋂EC))/(EC-(SC⋂EC))) with SC and EC the structural and effective connectivity matrixes, 
% and a ratio between 0 (equal amount of structural and effective connections) and ±∞ (zero structural or zero effective connections).

ratio = zeros(size(dataBase,2,1)); 
for subj=1:size(dataBase,2)

intersecting = dataBase(subj).dwi.SC_matrix - dataBase(subj).ccep.EC_matrix;

ratio(subj) = log10(sum(sum(intersecting == 1))/sum(sum(intersecting == -1))); %structural divided by effective
end

% change the order of the subject such that ecog are the first 5
j =1;
for i = [1,6,8,9,10,2,3,4,5,7,11,12,13] 
    ratio_streef(j,1) = ratio(i);
    j=1+j;
end
clear ratio intersecting j i subj 

V = figure('Renderer', 'painters', 'Position', [600 300 600 400]);
w = [1.02 1 0.98 1 1 ]; % make a XJitter such that the individual points do not overlap
scatter(w,ratio_streef(1:5),40,'o','k')
hold on
v = [0.74 0.7 0.7 0.7 0.7 0.6 0.7 0.67]; % make a XJitter such that the individual points do not overlap
scatter(v,ratio_streef(6:13),400,'.','k','LineWidth', 1.8,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
hold on 
yline(0)
ylim([-1.1 1.1])
xlim([-1 2])

saveas(V,'ratio','epsc') % save the figure for further processing with Adobe Illustrator

clear ratio_streef JI JI_alternative JI_expected JI_parker v w V

%% R warning
warning('run now R code in STReEF05_compare_networks_R and run section 1)')
%% SECTION 4: calculate the network topography and node proximity
%%
% this section was first implemented in the figures themselfs (not handy)
% change such that the network topography measures are also stored in the dataBase

% transform SOZ include such that it matches indexces of ch_include for
 % easy computation of the network measures
for subj = 1:size(dataBase,2)
ch = dataBase(subj).ccep.ch;
elec_include = dataBase(subj).dwi.elec_include;
soz = dataBase(subj).soz.soz;
ch_include = ch(elec_include);
soz_include = ch(soz);
soz_ind = NaN(size(soz_include,1),1);
for s = 1:size(soz_include,1)
    soz_in = soz_include(s);
    indx_soz = strcmpi(soz_in,ch_include); 
    if isempty(find(indx_soz,1))
        warning('soz channel is not included!') % toegevoegd 22/12/23 check nog waarom dit zo is
    else
    soz_ind(s) = find(indx_soz); 
    end
end
soz_ind = soz_ind(~isnan(soz_ind)); % toegevoegd 22/12/23
dataBase(subj).soz.soz_ind = soz_ind;
no_soz_ind = 1:size(ch_include,1);
no_soz_ind(soz_ind) = NaN;
no_soz_ind(isnan(no_soz_ind)) = [];
dataBase(subj).soz.no_soz_ind = no_soz_ind;
end

clear ch ch_include elec_include indx_soz s soz soz_in soz_include soz_ind subj no_soz_ind

disp('SOZ channels defined')
%% calculate network measures
for subj = 1:size(dataBase,2)

sub_label = erase(dataBase(subj).sub_label,'sub-');
SC = dataBase(subj).dwi.SC_matrix;
EC =  dataBase(subj).ccep.EC_matrix;  

ch = dataBase(subj).ccep.ch;
elec_include = dataBase(subj).dwi.elec_include;
ch_include = ch(elec_include);
soz_ind = dataBase(subj).soz.soz_ind;
no_soz_ind = dataBase(subj).soz.no_soz_ind;

% DEGREE
G_SC = graph(SC,ch_include); % graph structure
G_EC = graph(EC,ch_include);

degree_SC = degree(G_SC);
degree_EC = degree(G_EC);
dataBase(subj).network.degree_SC = degree_SC;
dataBase(subj).network.degree_EC = degree_EC;
dataBase(subj).network.degree_SC_soz = degree_SC(soz_ind);
dataBase(subj).network.degree_EC_soz = degree_EC(soz_ind);
dataBase(subj).network.degree_SC_nsoz = degree_SC(no_soz_ind);
dataBase(subj).network.degree_EC_nsoz = degree_EC(no_soz_ind);

% betweennes centrality
BC_SC = centrality(G_SC,'betweenness');
BC_EC = centrality(G_EC,'betweenness');
dataBase(subj).network.BC_SC = BC_SC;
dataBase(subj).network.BC_EC = BC_EC;
dataBase(subj).network.BC_SC_soz = BC_SC(soz_ind);
dataBase(subj).network.BC_EC_soz = BC_EC(soz_ind);
dataBase(subj).network.BC_SC_nsoz = BC_SC(no_soz_ind);
dataBase(subj).network.BC_EC_nsoz = BC_EC(no_soz_ind);
end
clear G_SC G_EC ch elec_include ch_include soz_ind no_soz_ind SC EC sub_label degree_SC degree_EC BC_SC BC_EC
%% node proximity
% calculate distance between two electrode contact coordinates
for subj = 1:size(dataBase,2)
elec_indx = dataBase(subj).dwi.elec_include;
sub_label = dataBase(subj).sub_label;
coordinates = [dataBase(subj).tb_electrodes.x(elec_indx) dataBase(subj).tb_electrodes.y(elec_indx) ...
                dataBase(subj).tb_electrodes.z(elec_indx)]; % the columns are the x,y,z coordinates, the rows the channels

% for some patients the variable tb_electrodes is stored in a different type of array. If so, change.
if isa(coordinates,'cell')
coordinates = str2double(coordinates);
end

distances = NaN(size(coordinates,1)); % distance in mm between two electrode contacts
for  cor = 1:size(coordinates,1)
for  cor2 = 1:size(coordinates,1)
D_dis = pdist2(coordinates(cor,:),coordinates(cor2,:),'euclidean'); % compute the euclidian distance between the electrode contact coordinates
if ~isempty(D_dis)
    distances(cor,cor2) = D_dis;
else 
    fprintf('empty coordinate in %d or %d',cor,cor2)
end
end
end
dataBase(subj).network.distances = distances; % save in dataBase
end

% compute the node proximity:  the node proximity per node was defined as the median distance between this node and all other nodes. 
for subj = 1:size(dataBase,2)
ch = dataBase(subj).ccep.ch;elec_include = dataBase(subj).dwi.elec_include;ch_include = ch(elec_include);
sub_label = dataBase(subj).sub_label;

node_prox = NaN(size(ch_include,1),1);

for ch = 1:size(ch_include,1)
node_prox(ch) = median(nonzeros(dataBase(subj).network.distances(:,ch)));
end
dataBase(subj).network.node_proximity = node_prox; % the node proximity per node was defined as the median distance between this node and all other nodes. 
end
clear ch_include ch coordinates cor cor2 D_dis distances elec_include elec_indx I_dis node_prox
%% prepare data for multilevel model
included =  NaN(size(cfg.sub_label,2),1);
for subj = 1:size(dataBase,2)
included(subj) = size(dataBase(subj).ccep.EC_matrix,1);
end
size_long = sum(included);
data_long = NaN(size_long,7); % matrix with for every channel per patient the specifications

i=1;
for subj = 1 %[1,2,3,4,6,7,8,9,10,11,12]% only patients where soz could be defined
sub_label = erase(dataBase(subj).sub_label,'sub-RESP');
SC = dataBase(subj).dwi.SC_matrix;
EC =  dataBase(subj).ccep.EC_matrix;  

ch = dataBase(subj).ccep.ch;
elec_include = dataBase(subj).dwi.elec_include;
ch_include = ch(elec_include);
sz = size(ch_include,1);
channels = 1:sz;
soz_ind = dataBase(subj).soz.soz_ind;
no_soz_ind = dataBase(subj).soz.no_soz_ind;

volume = dataBase(subj).dwi.volume_roi; % volume per electrode contact area
epi = NaN(sz,1);
epi(soz_ind) = 1; % if the channel is in the SOZ or not
epi(no_soz_ind) = 0;

data_long(i:i+sz-1,1)= str2double(sub_label)*ones(sz,1); % patient RESP numbers
data_long(i:i+sz-1,2)= dataBase(subj).network.degree_SC; % degree structural networks
data_long(i:i+sz-1,3)= dataBase(subj).network.degree_EC; % degree effective networks
data_long(i:i+sz-1,4)= dataBase(subj).network.node_proximity;
data_long(i:i+sz-1,5)= volume; % volume per electrode contact area
data_long(i:i+sz-1,6)= epi; % if the channel is in the SOZ or not
data_long(i:i+sz-1,7)= channels; % channel names

i=i+sz;
end

par_data = [" subj","degreeSC","degreeEC",'node_proximity','volume','epi','channels'];
data_l=array2table(data_long,'VariableNames',par_data);
dataPath = myDataPath.DWIMATLABpath; 
dir = [dataPath];
writetable(data_l,[dir 'data_long.csv']) % save the matrix with for every channel per patient the specifications

clear ch ch_include channels data_l data_long EC elec_include epi i included no_soz_ind par_data SC size_long soz_ind sz volume dir dataPath sub_label subj
%% R warning
warning('run now R code in STReEF05_compare_networks_R and run section 2)')

%% SECTION 5 save the important variables of the  file
% put all the database variables in one struct and save it as networkInfo.mat 
for subj = 1:size(dataBase,2)
    targetFolder = fullfile(myDataPath.DWIMATLABpath,dataBase(subj).sub_label,dataBase(subj).ses_label); % save the computed variables for further use in a designated folder 'dwi_matlab' 

    % Create the folder if it doesn't exist already.
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end
    
    start_filename = strfind(dataBase(subj).dwi.dataName,'/');
    stop_filename = strfind(dataBase(subj).dwi.dataName,'_ieeg');
    
    fileName = [dataBase(subj).dwi.dataName(start_filename(end)+1:stop_filename-1),'_networkInfo.mat']; % name it in a similair way as the CCEP files
    
    % insert variables
    network = dataBase(subj).network; % containing degree, BC, and node proximity 
    network.soz = dataBase(subj).soz; % containing soz and derivatives containing indexces of ch_include 
    network.SC_matrix = dataBase(subj).dwi.SC_matrix; % structural network matrix
    network.SC_matrix_values = dataBase(subj).dwi.SC_matrix_values; % structural network matrix with streamline density
    network.EC_matrix = dataBase(subj).ccep.EC_matrix; % structural network matrix
    network.EC_matrix_amplitude = dataBase(subj).ccep.EC_matrix_amplitude ; % effective network matrix with N1 amplitude
    network.EC_matrix_sample = dataBase(subj).ccep.EC_matrix_sample ; % effective network matrix with latency 

    save(fullfile(targetFolder,fileName), '-struct','network');
    fprintf('Saved network-struct in %s \n',fullfile(targetFolder,fileName))
end
clear fileName start_filename stop_filename targetFolder


 %% SECTION 6: functions to be implemented and check of the data
 %% exclude electrodes contact areas with too much overlap ( to be implemented in the analysis)

 for subj=1:size(dataBase,2)
 volume = dataBase(subj).dwi.volume_roi; %  volume per electrode contact area
  k=0;
 for i = 1:size(volume,1)
 if volume(i)<16 % if volume is under 16 mm3, exclude
     k = k+1;
 exclude(subj,k) = volume(i);
 exclude_ind(subj,k) = i;

 end
 end
 end

%% count the times a channel is stimulated
for subj=1:size(dataBase,2)
ch = dataBase(subj).ccep.ch;
elec_include = dataBase(subj).dwi.elec_include;
ch_include = ch(elec_include);
cc_stimchans = dataBase(subj).ccep.cc_stimchans;
stimpair = NaN(size(ch_include,1),1);

for i = 1:size(ch_include,1)
    chan = ch_include(i);
    loc = strcmpi(chan,cc_stimchans);
    [col,~] = find(loc); 
    stimpair(i) = size(col,1);
end
dataBase(subj).stimcount = stimpair;
end
clear stimpair col loc chan i ch_include cc_stimchans ch elec_include subj
% calculate the unique possibilities of the amount of times a channel is stimulated
for subj=1:size(dataBase,2)
    stimsorts=unique(dataBase(subj).stimcount);
    dataBase(subj).stim_sort = stimsorts;
end
clear stimsorts 


