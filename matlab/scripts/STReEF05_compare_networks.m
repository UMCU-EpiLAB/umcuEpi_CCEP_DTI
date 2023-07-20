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

%% load SOZ information 
% Run this only once and save (information is loaded in next section)


subj = 1; % run this part for every subject
ch = dataBase(subj).ccep.ch;
elec_include = dataBase(subj).dwi.elec_include;
soz_include = false(size(elec_include));

% insert the SOZ as defined by 'bijlage 3' by hand ((we just have to insert this in the BIDS now I know how to do that)
x1 = input(sprintf('Insert SOZ channels of subj %s:',dataBase(subj).sub_label),'s');
x = strsplit(x1,{', ',','});
for i = 1:size(x,2)
x2 = x(i);
loc = strcmpi(x2,ch);
[ch_ind,~] = find(loc); % gaat om de row
if ~isempty(ch_ind)
if elec_include(ch_ind)
    soz_include(ch_ind) = 1;
else disp(sprintf('Channel %s not included',string(cellstr(x2)))) % warn if there is a channel not included
end
else disp(sprintf('Channel %s not found',string(cellstr(x2)))) % warn if there is a channel not found
end
end
dataBase(subj).soz = soz_include;

% save SOZ
for subj = 1:size(dataBase,2)
    targetFolder = fullfile(myDataPath.DWIMATLABpath,dataBase(subj).sub_label,dataBase(subj).ses_label);

    % Create the folder if it doesn't exist already.
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end
    
    start_filename = strfind(dataBase(subj).dwi.dataName,'/');
    stop_filename = strfind(dataBase(subj).dwi.dataName,'_ieeg');
    
    fileName = [dataBase(subj).dwi.dataName(start_filename(end)+1:stop_filename-1),'_sozInfo.mat'];
            
    soz = struct();
    soz.dataName = dataBase(subj).dwi.dataName;
    soz.soz_include = dataBase(subj).soz;
   
    save(fullfile(targetFolder,fileName), '-struct','soz');
    fprintf('Saved soz in %s \n',fullfile(targetFolder,fileName))
end

%% load processed electrode contact area data (see STReEF02_coreg_roidef_mrtrix) and visual scored CCEP data (see umcuEpi_CCEP_DTI (main)/matlab/scripts/ccepDTI03_detectCCEP.m for scoring CCEP data)
% merge runs: if size(cfg.run_label{subj},2) > 1
dataPath = myDataPath.CCEP;
dataPath2 = myDataPath.DWIMATLABpath;
dataBase = struct([]);
for subj = 1:size(cfg.sub_label,2)
    sub_label = ['sub-' cfg.sub_label{subj}];
    ses_label = cfg.ses_label{subj};
    task_label = cfg.task_label{subj};
    run_label_dwi = cfg.run_label{subj}{1}; % er is maar een dwi per subj, zo is het opgeslagen beetje onhandig
    
    % load dwiInfo + sozInfo which contains information about the included electrode contact coordinates and the electrode contact areas
    DWI = dir(fullfile(dataPath2,sub_label,ses_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label_dwi,'_dwiInfo.mat'])); % made by running STReEF02_coreg_roidef_matlab
    SOZ = dir(fullfile(dataPath2,sub_label,ses_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label_dwi,'_sozInfo.mat'])); % load SOZ if you runned the 'load SOZ information' earlier
    if size(DWI,1) == 0
            error('%s does not exist',fullfile(dataPath2,sub_label,ses_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label_dwi,'_dwiInfo.mat']));
    end 
    if size(SOZ,1) == 0
            error('%s does not exist',fullfile(dataPath2,sub_label,ses_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label_dwi,'_sozInfo.mat']));
    end 
    dataNameDWI = fullfile(DWI(1).folder, DWI(1).name);
    dwi = load(dataNameDWI);
    dataNameSOZ = fullfile(SOZ(1).folder, SOZ(1).name);
    soz = load(dataNameSOZ);

    dataBase(subj).ses_label = ses_label;
    dataBase(subj).sub_label = sub_label;
    dataBase(subj).dwi = dwi;
    dataBase(subj).soz = soz.soz_include;

    % load detected and visual scored ccep data (made by running umcuEpi_CCEP_DTI (main)/matlab/scripts/ccepDTI03_detectCCEP.m)
    for run = 1:size(cfg.run_label{subj},2)
    run_label = cfg.run_label{subj}{run}; 
 
    DRE = dir(fullfile(dataPath,sub_label,ses_label,run_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label,'_N1sChecked.mat'])); % checked data
    if size(DRE,1) == 0
            error('%s does not exist',fullfile(dataPath,sub_label,ses_label,run_label,...
            [sub_label, '_', ses_label,'_',task_label,'_',run_label,'_N1sChecked.mat']));
    end 
    dataNameRE = fullfile(DRE(1).folder, DRE(1).name);
    ccep = load(dataNameRE);

    dataBase(subj).ccep(run).run_label = run_label;
    dataBase(subj).ccep(run).ccep = ccep;
    end
end
clear dataNameRE dataPath run_label ses_label sub_label task_label DRE run_label_dwi run subj
disp('visual checked data loaded')
clear ccep dataNameDWI dataNameSOZ dataPath2 dwi DWI soz SOZ
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
clear scored run subj stimsets stimchannels  n1_peak_sample n1_peak_amplitude
%% solve mistake in stimchans of patient 10 (10/13) (we just have to change this in the BIDS now I know how to do that)
 dataBase(10).ccep.cc_stimchans(48,2)={'F07'};
 %% SECTION 2: construct the effective and structural networks
%% construct an effective network matrix
% Connections were drawn from both electrode contacts in a stimulus pair to the electrode contacts in which an ER was detected.  

for subj=1:size(dataBase,2)
    checked = dataBase(subj).ccep.checked; % matrix with zeros and ones for detected (and visual checked) ERs
    n1_peak_amplitude = dataBase(subj).ccep.n1_peak_amplitude; % matrix with zeros and values of n1 peak amplitude
    n1_peak_sample = dataBase(subj).ccep.n1_peak_sample; % matrix with zeros and peak sample value of n1 peak amplitude
    n1_peak_sample(isinf(n1_peak_sample)) = NaN; % make the inf NaN (Inf are in the matrices due to pre-allocation in visualRating_ccep)
    n1_peak_amplitude(isinf(n1_peak_amplitude)) = NaN; % make the inf NaN

    ch = dataBase(subj).ccep.ch; % all channel names
    elec_include = dataBase(subj).dwi.elec_include; % logical vector including all channels and their inclusion yes (1) or no (0) for easy computational matters
    cc_stimchans = dataBase(subj).ccep.cc_stimchans; % channel names of stimulated electrode pairs
    cc_stimsets = dataBase(subj).ccep.cc_stimsets;  % channel number of stimulated electrode pairs
    ch_include = ch(elec_include); % channel names of the included channels

    EC_matrix_ccep = zeros(size(find(elec_include),1)); % effective matrix with zeros and ones with on the rows and columns the channels
    EC_matrix_amplitude = NaN(size(find(elec_include),1)); % effective matrix with n1 peak amplitude
    EC_matrix_sample = NaN(size(find(elec_include),1)); % effective matrix with n1 peak sample value

    k=0;
    % run to every stimulated pair to place the detected ERs in the right spot in the matrix with on the rows and columns each included channel one time
    for i = 1:size(find(elec_include),1)
    chan = ch_include(i); % channel name
    loc = strcmpi(chan,cc_stimchans); % stim pairs were this channel was included in
    [col,~] = find(loc); 
    if size(col,1)>1 % if channel was included in two stimulation pairs, combine these results
        col_data = sum(checked(elec_include,col),2); % define a connection only if there is a ER evoked in after stimulating both stimulation pairs NaN+NaN=NaN, Inf+Inf = Inf Nan+Inf = NaN 1+Inf = NaN 1+NaN = NaN 1+1 = 2
        [col2,~] = find(checked(elec_include,col)==1); % to all 1+Inf = NaN 1+NaN = NaN 1+1 = 2 assign 1
        col_data(col2,1)=1;
        col_data_amplitude = mean(n1_peak_amplitude(elec_include,col),2,'omitnan');
        col_data_sample = mean(n1_peak_sample(elec_include,col),2,'omitnan');
    else
        col_data = checked(elec_include,col); % extract the responses from the included channels after stimulating the selected channel
        col_data_amplitude = n1_peak_amplitude(elec_include,col); % extract the n1 peak amplitude
        col_data_sample = n1_peak_sample(elec_include,col); % extract the n1 peak sample value
    end
    
    if ~isempty(col_data)
    EC_matrix_ccep(:,i) = col_data; % insert the response in the effective matrix
    EC_matrix_amplitude(col_data==1,i) = col_data_amplitude(col_data==1); % insert the n1 peak amplitude for only the detected (and visually checked) ERs
    EC_matrix_sample(col_data==1,i) = col_data_sample(col_data==1); % insert the n1 peak sample value for only the detected (and visually checked) ERs
    end

    % deal with included channels that are not stimulated on by making the response after stimulating all other channels the virtual response after stimulating this channel
    if isempty(col_data)
    k = k + 1;
    no_stim_col(k) = i;
    end
    end

    if exist('no_stim_col') 
    for j = 1:size(no_stim_col,2)
    chan = no_stim_col(j);
    col_data = EC_matrix_ccep(chan,:); % for channels not included in a stim pair, make de ERs evoked after stimulating all other channels the response in the effective matrix (since we are making an symmetric matrix later on, we can do this)
    EC_matrix_ccep(:,chan) = col_data;
    EC_matrix_amplitude(:,chan) = EC_matrix_amplitude(chan,:);
    EC_matrix_sample(:,chan) = EC_matrix_sample(chan,:);
    end
    end

    EC_matrix_ccep(isnan(EC_matrix_ccep))=0; %make all the NaN zero's
    EC_matrix_ccep(EC_matrix_ccep==Inf)=0; %make al the Inf zero's
    dataBase(subj).ccep.EC_matrix = EC_matrix_ccep; % store the effective matrices in the dataBase for later use
    dataBase(subj).ccep.EC_matrix_amplitude = EC_matrix_amplitude;
    dataBase(subj).ccep.EC_matrix_sample = EC_matrix_sample; 
    clear('no_stim_col')
end

disp('EC_matrix made')

clear EC_matrix_amplitude EC_matrix_sample n1_peak_amplitude n1_peak_sample col_data_amplitude col_data_sample data run_label cc_stimchans cc_stimsets ccep ch ch_include chan checked col_data dwi DWI elec_include i j k loc col col2 EC_matrix_ccep subj

%% make the effective network matrix symmetric
% Effective networks were made symmetrical by considering all ERs as bi-directional, to be able to compare to the non-directional structural networks.
for subj=1:size(dataBase,2)
    EC_matrix = dataBase(subj).ccep.EC_matrix;
    EC_matrix_amplitude = dataBase(subj).ccep.EC_matrix_amplitude;
    EC_matrix_sample  = dataBase(subj).ccep.EC_matrix_sample;

    EC_matrix_sym = NaN(size(EC_matrix,1));
    EC_matrix_sym_amplitude = NaN(size(EC_matrix,1));    
    EC_matrix_sym_sample = NaN(size(EC_matrix,1));  

    for i = 1:size(EC_matrix,1)
    connections = EC_matrix(:,i)|EC_matrix(i,:)'; % consider an connection between channel i and channel j if there is a ER evoked in either direction (so after stimulating channel i or channel j) 
    EC_matrix_sym(:,i) = connections;
    EC_matrix_sym(i,:) = connections;
    
    both_amplitude = [EC_matrix_amplitude(connections,i),EC_matrix_amplitude(i,connections)'];
    EC_matrix_sym_amplitude(connections,i) = mean(abs(both_amplitude),2,'omitnan'); % consider the mean amplitude (omiting NaNs) of the detected ERs in both directions
    EC_matrix_sym_amplitude(i,connections) =  mean(abs(both_amplitude),2,'omitnan');

    both_sample = [EC_matrix_sample(connections,i),EC_matrix_sample(i,connections)']; % consider the mean sample value (omiting NaNs) of the detected ERs in both directions
    EC_matrix_sym_sample(connections,i) = mean(abs(both_sample),2,'omitnan');
    EC_matrix_sym_sample(i,connections) =  mean(abs(both_sample),2,'omitnan');
    end
    dataBase(subj).ccep.EC_matrix = EC_matrix_sym; 
    dataBase(subj).ccep.EC_matrix_amplitude = EC_matrix_sym_amplitude;
    dataBase(subj).ccep.EC_matrix_sample = EC_matrix_sym_sample;     
end
clear EC_matrix EC_matrix_sym_sample EC_matrix_sym_amplitude EC_matrix_sym EC_matrix_amplitude EC_matrix_sample EC_matrix_sym i subj connections both_sample both_amplitude
disp('symmetric matrix made')


%% construct an structural network matrix
% load in the connectome matrix (preliminary structural network) made with STReEF04_fiber_tractography
for subj=1:size(dataBase,2)
dataPath = myDataPath.FTpath; 
sub_label = erase(dataBase(subj).sub_label,'sub-RESP');
dir1 = [dataPath dataBase(subj).sub_label];

% define tractography parameters to open the right file with the correct parameters
treshold = 0.1; % a structural connection was formed when the streamline density exceeded a threshold of 0.1.
cutoff = 0.15;
seed_voxel = 6000; 
angle = 70; 
step = 1;  
minlength = 4; 
maxlength = 400;

dataName = sprintf('/connectome_backtrack_%s_%d_%d_%d_%d_%d_%g_%d.csv',...
    sub_label,seed_voxel,step,angle,minlength,maxlength,cutoff,angle);
SC_matrix_dwi = load([dir1 dataName]);

% check if the sizes of the effective and prelimary structural network are the same
if size(SC_matrix_dwi,2) ~= size(dataBase(subj).ccep.EC_matrix,2)
disp(subj)
end

dataBase(subj).dwi.SC_matrix_values = SC_matrix_dwi; 

SC_matrix_dwi(SC_matrix_dwi<treshold)=0; % a structural connection was only formed when the streamline density exceeded a threshold of 0.1.
SC_matrix_dwi(SC_matrix_dwi>=treshold)=1; 

dataBase(subj).dwi.SC_matrix = SC_matrix_dwi; 
end
disp('SC_matrix made')
 clear angle cutoff maxlength minlength SC_matrix_dwi seed_voxel step treshold cor_size dir1 dataPath dataName subj sub_label

 %% define the SOZ channels correctly
 % transform SOZ include such that it matches indexces of ch_include
for subj = 1:size(dataBase,2)
ch = dataBase(subj).ccep.ch;
elec_include = dataBase(subj).dwi.elec_include;
soz = dataBase(subj).soz;
ch_include = ch(elec_include);
soz_include = ch(soz);
soz_ind = NaN(size(soz_include,1),1);
for s = 1:size(soz_include,1)
    soz_in = soz_include(s);
    indx_soz = strcmpi(soz_in,ch_include); 
    soz_ind(s) = find(indx_soz); 
end
dataBase(subj).soz_ind = soz_ind;
no_soz_ind = 1:size(ch_include,1);
no_soz_ind(soz_ind) = NaN;
no_soz_ind(isnan(no_soz_ind)) = [];
dataBase(subj).no_soz_ind = no_soz_ind;
end

clear ch ch_include elec_include indx_soz s soz soz_in soz_include soz_ind subj no_soz_ind

% find the indeces from the vector with all channel names for the included and soz channels for easier computations
for subj = 1:size(dataBase,2)
ch = dataBase(subj).ccep.ch;
elec_include = dataBase(subj).dwi.elec_include;
ch_include = ch(elec_include);
soz = dataBase(subj).soz;
soz_include = ch(soz); 

elec_ind = NaN(size(ch_include,1),1);
soz_ind_ch = NaN(size(ch(soz),1),1);
for s = 1:size(ch_include,1)
    elec_in = ch_include(s);
    indx = strcmpi(elec_in,ch); 
    elec_ind(s) = find(indx); 
end
for s = 1:size(soz_include,1)
    soz_in = soz_include(s);
    indx = strcmpi(soz_in,ch);
    soz_ind_ch(s) = find(indx);
end
dataBase(subj).soz_ind_ch = soz_ind_ch; % indeces of the soz channels when counting with all channels
dataBase(subj).elec_ind = elec_ind; % indeces of the included channels when counting with all channels
end
clear ch ch_include elec_include indx s elec_in elec_ind subj soz_ind_ch soz soz_in soz_include 

%% SECTION 3: caculate the inter-modal similarity and prepare for the Jaccard Index calculation in R

JI = zeros(size(dataBase,2,1));
for subj=1:size(dataBase,2)
sub_label = erase(dataBase(subj).sub_label,'sub-RESP');
dir_SC = ['/Fridge/users/susanne/derivatives/dwi_matlab/' dataBase(subj).sub_label];
a = dataBase(subj).dwi.SC_matrix;
b =  dataBase(subj).ccep.EC_matrix;  
save([dir_SC '/SC_matrix_dwi'],'a','-ascii') % save the structural network
save([dir_SC '/SC_matrix_ccep'],'b','-ascii') % save the effective network
JI(subj) = sum(sum(a & b))/sum(sum(a | b)); % calculate the Jaccard Index to be able to compare it with the Jaccard Index calculated with R
end
sub_labels = cfg.sub_label;
grid = [1,6,8,9,10];
JI_grid = table(sub_labels(grid)',JI(grid));
seeg = [2,3,4,5,7,11,12,13];
JI_seeg = table(sub_labels(seeg)',JI(seeg));

median_seeg = median(JI_seeg.Var2)
IQR_seeg = iqr(JI_seeg.Var2)
median_grid = median(JI_grid.Var2)
IQR_grid = iqr(JI_grid.Var2)
prctile(JI_grid.Var2,[25 75])
prctile(JI_seeg.Var2,[25 75])

result_JI = table( [round(JI,2);median(JI);prctile(JI,[25 75])']); % present the statisctics of the Jaccard Index in a table

clear grid seeg sub_label sub_labels subj dir_SC a b IQR_grid IQR_seeg median_grid median_seeg %JI
%% R warning
warning('run now R code in STReEF05_compare_networks_R and run section 1)')
%% SECTION 4: calculate the network topography
%%
% this section was first implemented in the figures themselfs (not handy)
% change such that the network topography measures are also stored in the dataBase
for subj = 1:size(dataBase,2)

sub_label = erase(dataBase(subj).sub_label,'sub-');
SC = dataBase(subj).dwi.SC_matrix;
EC =  dataBase(subj).ccep.EC_matrix;  

ch = dataBase(subj).ccep.ch;
elec_include = dataBase(subj).dwi.elec_include;
ch_include = ch(elec_include);
soz_ind = dataBase(subj).soz_ind;
no_soz_ind = dataBase(subj).no_soz_ind;

% DEGREE
G_SC = graph(SC,ch_include); % graph structure
G_EC = graph(EC,ch_include);

degree_SC = degree(G_SC);
degree_EC = degree(G_EC);
degree_SC_soz = degree_SC(soz_ind);
degree_EC_soz = degree_EC(soz_ind);
degree_SC_nsoz = degree_SC(no_soz_ind);
degree_EC_nsoz = degree_EC(no_soz_ind);

% betweennes centrality
BC_SC = centrality(G_SC,'betweenness');
BC_EC = centrality(G_EC,'betweenness');
BC_SC_soz = BC_SC(soz_ind);
BC_EC_soz = BC_EC(soz_ind);
BC_SC_nsoz = BC_SC(no_soz_ind);
BC_EC_nsoz = BC_EC(no_soz_ind);
end
%% prepare data for multilevel model

data_long = NaN(size(d_all,1),11); % matrix with for every channel per patient the specifications

i=1;
for subj = [1,2,3,4,6,7,8,9,10,11,12]% only seizure free patients
SC = dataBase(subj).dwi.SC_matrix;
EC =  dataBase(subj).ccep.EC_matrix;  
sub_label = erase(dataBase(subj).sub_label,'sub-RESP');
ch = dataBase(subj).ccep.ch;
elec_include = dataBase(subj).dwi.elec_include;
ch_include = ch(elec_include);
sz = size(ch_include,1);
channels = 1:sz;
soz_ind = dataBase(subj).soz_ind;
no_soz_ind = dataBase(subj).no_soz_ind;
distance = dataBase(subj).dwi.distance_ch;
center_distance = distance - median(distance); % node proximity
volume = dataBase(subj).dwi.volume_roi; % volume per electrode contact area
epi = NaN(sz,1);
epi(soz_ind) = 1; % if the channel is in the SOZ or not
epi(no_soz_ind) = 0;

% DEGREE
G_SC = graph(SC,ch_include); % graph structure
G_EC = graph(EC,ch_include);

degree_SC = degree(G_SC);
degree_EC = degree(G_EC);
center_SC = degree_SC - median(degree_SC);
center_EC = degree_EC - median(degree_EC);
nor_SC = degree_SC./sz;
nor_EC = degree_EC./sz;

data_long(i:i+sz-1,1)= str2double(sub_label)*ones(sz,1); % patient RESP numbers
data_long(i:i+sz-1,2)= degree_SC; % degree structural networks
data_long(i:i+sz-1,3)= degree_EC; % degree effective networks
data_long(i:i+sz-1,4)= center_SC; % centered degree (not used in the end)
data_long(i:i+sz-1,5)= center_EC;
data_long(i:i+sz-1,6)= nor_SC; % normalized degree with the nr of channels (not used in the end)
data_long(i:i+sz-1,7)= nor_EC;
data_long(i:i+sz-1,8)= distance;
data_long(i:i+sz-1,9)= center_distance; % node proximity
data_long(i:i+sz-1,10)= volume; % volume per electrode contact area
data_long(i:i+sz-1,11)= epi; % if the channel is in the SOZ or not
data_long(i:i+sz-1,12)= channels; % channel names

i=i+sz;
end
par_data = [" subj","degreeSC","degreeEC",'centerSC','centerEC','nor_SC','nor_EC','distance','centerdistance','volume','epi','channels'];
data_l=array2table(data_long,'VariableNames',par_data);
writetable(data_l,'data_long.csv') % save the matrix with for every channel per patient the specifications
%% R warning
warning('run now R code in STReEF05_compare_networks_R and run section 2)')

%% SECTION 5: Node proximity
% calculate distance between two electrode contact coordinates (this section can be shorter)
% this section was first implemented in the figures themselfs (not handy)
% change such that the node proximity are also stored in the dataBas
for subj = 1:size(dataBase,2)
elec_indx = dataBase(subj).metadata.elec_include;
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
[D_dis, I_dis] = pdist2(coordinates(cor),coordinates(cor),'euclidean','Largest',1); % compute the euclidian distance between the electrode contact coordinates
if ~isempty(D_dis)
    distances(cor,cor2) = max(D_dis);
else 
    fprintf('empty roi volume van %d or %d',cor,cor2)
end
end
end
dataBase(subj).metadata.distances_ch = distances; % save in dataBase

dir2 = [myDataPath.coreg_ROIpath sprintf('%s/',sub_label)];
save([dir2 'distances_ch.txt'],'distances','-ascii') % save in the designated co-registration and roi definition folder ('coreg_ROI') 
end
% compute the node proximity:  the node proximity per node was defined as the median distance between this node and all other nodes. 
for subj = 1:size(dataBase,2)
ch = dataBase(subj).ccep.ch;elec_include = dataBase(subj).dwi.elec_include;ch_include = ch(elec_include);
sub_label = dataBase(subj).sub_label;

distances_median = NaN(size(ch_include,1),1);

dir2 = [myDataPath.coreg_ROIpath sprintf('%s/',sub_label)];
distance = load([dir2 'distances_ch.txt']);

for ch = 1:size(ch_include,1)
distances_median(ch) = median(nonzeros(distance(:,ch)));
end
dataBase(subj).dwi.distance_ch = distances_median; % the node proximity per node was defined as the median distance between this node and all other nodes. 
% node proximity = center_distance = distance - median(distance);
end

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


