% STReEF06_figures
% make figures to present the inter-modal similarity and network topography comparison between structural (derived from DWI) and effective (derived from SPES) networks

% author: Susanne Jelsma & Dorien van Blooijs
% date: June 2022

%% SECTION 1: prepare data 

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

%% load processed network data(see STReEF05_compare_networks) electrode contact area data (see STReEF02_coreg_roidef_mrtrix) and visual scored CCEP data (see umcuEpi_CCEP_DTI (main)/matlab/scripts/ccepDTI03_detectCCEP.m for scoring CCEP data)

dataBase = load_derivatives_data(myDataPath,cfg);

disp('network data loaded')

%% merge runs 
% Be aware! Code with potential of causing errors due to copying of information, so check for new types of patients or new code added!
% implement in function
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

%% SECTION 2  
% visualize the inter-modal similarity with connectivity matrices and visualize the network topography per patient

%% visualize the inter-modal similarity with connectivity matrices

for subj = 1:size(dataBase,2)
V = visual_networks(dataBase(subj).network.SC_matrix,dataBase(subj).network.EC_matrix,dataBase(subj).sub_label) % function to plot the connectivity matrices
V.WindowState = 'maximized';
print('-vector','-depsc',V,sprintf('visual_networks_symetric%s',dataBase(subj).sub_label))% save the figure for further processing with Adobe Illustrator
end

%%  visualize the network topography for all patients degree
i=0;j=0; 
for subj =[1,6,8,9,10]; %[2,3,4,5,7,11] %seeg1 %[12,13] seeg2 %[1,6,8,9,10]grid % 
i=i+1

V = visual_topography(dataBase(subj).network.degree_EC, dataBase(subj).network.degree_SC, dataBase(subj).network.degree_EC_soz, dataBase(subj).network.degree_EC_nsoz,dataBase(subj).network.degree_SC_soz, dataBase(subj).network.degree_SC_nsoz,subj,i)
end
saveas(V,'correlation degree grid','epsc') % save the figure for further processing with Adobe Illustrator

%%  visualize the network topography for all patients betweenness centrality
i=0;j=0; 
for subj = [1,6,8,9,10]; %[2,3,4,5,7,11] %seeg1 %[12,13] seeg2 %[1,6,8,9,10]grid %
i=i+1

V = visual_topography(dataBase(subj).network.BC_EC, dataBase(subj).network.BC_SC, dataBase(subj).network.BC_EC_soz, dataBase(subj).network.BC_EC_nsoz,dataBase(subj).network.BC_SC_soz, dataBase(subj).network.BC_SC_nsoz,subj,i)
end
saveas(V,'correlation BC seeg2','epsc') % save the figure for further processing with Adobe Illustrator

%%  visualize the network topography for all patients node proximity
i=0;j=0; 
for subj = [1,6,8,9,10] %1:size(dataBase,2)  [5,7,11,12]13[1,6,8,9][10,2,3,4]
i=i+1
np_soz = dataBase(subj).network.node_proximity(dataBase(subj).network.soz.soz_ind) ;
np_nsoz = dataBase(subj).network.node_proximity(dataBase(subj).network.soz.no_soz_ind) ; 
V = visual_topography_np(dataBase(subj).network.node_proximity, dataBase(subj).network.degree_EC, np_soz, np_nsoz ,dataBase(subj).network.degree_EC_soz, dataBase(subj).network.degree_EC_nsoz,subj,i)
end
%saveas(V,'correlation np effective grid','epsc') % save the figure for further processing with Adobe Illustrator

%%  visualize the network topography for all patients volume of electrode contact areas (VEA)
i=0;j=0; 
for subj = [12,13] %seeg2 %[1,6,8,9,10]grid [2,3,4,5,7,11] %seeg1 %
i=i+1
np_soz = dataBase(subj).dwi.volume_roi(dataBase(subj).network.soz.soz_ind) ;
np_nsoz = dataBase(subj).dwi.volume_roi(dataBase(subj).network.soz.no_soz_ind) ; 
V = visual_topography_np(dataBase(subj).dwi.volume_roi, dataBase(subj).network.degree_SC, np_soz, np_nsoz ,dataBase(subj).network.degree_SC_soz, dataBase(subj).network.degree_SC_nsoz,subj,i)
end
saveas(V,'correlation VEA structural seeg2','epsc') % save the figure for further processing with Adobe Illustrator

%% Histogram of the volume of electrode contact areas
i=0;j=0; 
for subj = [2,3,4,5,7,11] %[2,3,4,5,7,11] %seeg1 %[12,13] seeg2 %[1,6,8,9,10]grid
i=i+1;
set(gcf,'renderer','Painters')
V = figure(1);
V.WindowState = 'maximized';
 
subplot(2,3,i)

vea_soz = dataBase(subj).dwi.volume_roi(dataBase(subj).network.soz.soz_ind) ;
vea_nsoz = dataBase(subj).dwi.volume_roi(dataBase(subj).network.soz.no_soz_ind) ; 
histogram(vea_nsoz,'BinWidth',16,'BinLimits',[0,64],'FaceColor','k','FaceAlpha',1,"LineStyle","-") % plot volume electrode contact areas
axis square
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
hold on
if ~isempty(vea_soz)
histogram(vea_soz,'BinWidth',16,'BinLimits',[0,64],'FaceColor',[204/250 37/250 41/250],'FaceAlpha',1,'EdgeColor','none') % plot volume electrode contact areas in the soz
end
histogram(vea_nsoz,'BinWidth',16,'BinLimits',[0,64],'EdgeColor','k','DisplayStyle','stairs') % plot line histogram of non soz for when soz counts overlap the non soz counts

xlim([0,64])

set(gca, 'XTick',[0,16,32,48,64])
set(gca, 'XTickLabel',[0,16,32,48,64])
set(gca, 'YTick',[0,max(yt)/2,max(yt)])     
set(gca, 'YTickLabel',[0,max(yt)/2,max(yt)])
set(gca, 'FontSize',16) 
end
saveas(V,'histogram VEA seeg1','epsc') % save the figure for further processing with Adobe Illustrator

clear vea_nsoz w xt yt V vea_soz i j subj
