% STReEF02_coreg_roidef_mrtrix
% coregistration and determination of electrode contact coordinates for region-of-interest definition

% author: Susanne Jelsma & Dorien van Blooijs
% date: February 2022

% Be aware! This file has a twin written in mrtrix code and to execute the code correctly, the sections of these files must be runned alternately. 
% So first section 1 of the file: ' STReEF02_coreg_roidef_mrtrix' and than section 1 of this file so on.. 

% load sEEG/ECOG data, split into stimulation trials, select electrode channels, co-register MRI and DWI data and calculate the transformation matrix,
% transform electrode contact coordinates, define region-of-interests(ROI)/electrode contact areas, and save the important variables

%% SECTION 1: prepare sEEG/ECOG data (load sEEG/ECOG data, split into stimulation trials), select electrode channels, co-register MRI and DWI data with SPM, and calculate the transformation matrix  

%% open spm now in the command window if you want to check the co-registration later
%% set paths
% set umcuEpi_DTI/matlab in your directory and run this section

clc
clear
cfg.folderinput = 'chronic_ECoG'; % from which folder would you like to load ECoGs?
myDataPath = setLocalDataPath(cfg);

%% patient characteristics
%search for the SPES patients who have a DWI scan

files = dir(myDataPath.DWIpath); % path to the folder with the processed DWI's
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

%% load ECoGs 

dataBase = load_ECoGdata(myDataPath,cfg);

disp('All ECoGs are loaded')

%% preprocessing CCEP in ECoG
%clear cfg (removed this to remain some patient info (run_label etc) in the cfg we use
%later on)

% preprocessing step
cfg.dir = 'no'; % if you want to take negative/positive stimulation into account
cfg.amp = 'no'; % if you want to take stimulation current into account

% select epochs and average
cfg.epoch_length = 4; % in seconds, -2:2
cfg.epoch_prestim = 2;

dataBase = preprocess_ECoG(dataBase, cfg);

disp('All ECoGs are preprocessed')

%% rereference data
% I did not include the rereferencing step in the data processing (to shorten the run time) because in this file we are only looking for the electrode positions and
% not displaying the ecog data. I did the preprocessing step because I programmed the code beneath with the format of the dataBase struct after the
% preprocess_ECoG function (especially need the columns dataBase.metadata.{cc_stimchans and cc_stimsets}).
%                  
% Be aware! The rerefCCEPdata function is not necessary if looking at the electrode positions only but you need to use it when
% displaying ecog data as well! (merge also more things than in the next loop)
%% merge runs
% Be aware! Code with potential of causing errors due to copying of information, so check for new types of patients!

for subj= 1:size(dataBase,2)
    if size(cfg.run_label{subj},2) > 1
       dataBase(subj).metadata_runs = dataBase(subj).metadata; % copy everything to metadata_runs
       data = dataBase(subj).metadata(1); % start with the basic info from the first run
       dataBase(subj).metadata = data;
       % append the vectors cc_stimsets (containing channel number of stimulated electrode pairs) and cc_stimchans (containing channel names of stimulated electrode pairs) from each run into one large vector
            for run = 2:size(cfg.run_label{subj},2)
                stimsets = dataBase(subj).metadata_runs(run).cc_stimsets;
                stimchannels = dataBase(subj).metadata_runs(run).cc_stimchans;

                dataBase(subj).metadata.cc_stimsets = [dataBase(subj).metadata.cc_stimsets;stimsets];
                dataBase(subj).metadata.cc_stimchans = [dataBase(subj).metadata.cc_stimchans;stimchannels];
            end
        fprintf('...Runs subject %s has been merged ... \n',dataBase(subj).sub_label)
    end
end

disp('runs merged')
clear data run subj stimsets stimchannels stimchannels


%% select and check the electrode channels

for  subj=1:size(dataBase,2)
run=1;
% if sEEG, use only the gray matter channels, hippocampus, amygdala, lesion, and gliosis
if any(contains(fieldnames(dataBase(subj).tb_electrodes),'graymatter')) % select sEEG patients
    idx_screw = strcmpi(dataBase(subj).tb_electrodes.screw,'yes'); % remove screw electrodes
    idx_csf = strcmpi(dataBase(subj).tb_electrodes.csf,'yes'); % remove csf electrodes
    idx_whitematter = strcmpi(dataBase(subj).tb_electrodes.whitematter,'yes') & strcmpi(dataBase(subj).tb_electrodes.graymatter,'no'); % remove white matter channels (not the bordeline gray/white matter)
    idx_bad = strcmpi(dataBase(subj).metadata(run).tb_channels.status,'bad'); % remove bad channels from analysis because you cannot make extract ERs from it, so cannot compare it to the structural networks
    idx_all = sum([idx_screw, idx_csf, idx_whitematter, idx_bad],2);

    elec_include = cell(size(dataBase(subj).tb_electrodes,1),1);
    [elec_include{idx_all>0}] = deal('no'); % no in elec_include vector if it is a screw, csf, whitematter, or bad channel
    [elec_include{idx_all==0}] = deal('yes'); % yes in elec_include vector if it is something else (options: graymatter, hippocampus, amygdala, lesion, or gliosis)
elseif any(contains(dataBase(subj).tb_electrodes.group,'grid')) % select grid patients
    idx_silicon = strcmpi(dataBase(subj).tb_electrodes.silicon,'yes');% remove electrodes who are laying on other electrodes (silicon)
    idx_bad = strcmpi(dataBase(subj).metadata(run).tb_channels.status,'bad'); % remove bad channels from analysis because you cannot make extract ERs from it, so cannot compare it to the structural networks
    idx_all = sum([idx_silicon, idx_bad],2);

    bad = length(find(idx_bad==1))-length(find(idx_silicon==1));
    warning('%s extra bad channels %g',dataBase(subj).sub_label,bad) % warning because this was wrongly not included in earlier version (only 1 pt, 4 channels)

    elec_include = cell(size(dataBase(subj).tb_electrodes,1),1);
    [elec_include{idx_all>0}] = deal('no');
    [elec_include{idx_all==0}] = deal('yes');
else
    warning('%s is not a sEEG or grid',dataBase(subj).sub_label)
    elec_include = cell(size(dataBase(subj).tb_electrodes,1),1);
    [elec_include{1:end}] = deal('yes');
end

% save in the dataBase for further use
chan_include = dataBase(subj).metadata(run).ch(strcmpi(elec_include,'yes'));
dataBase(subj).metadata.ch_include = chan_include; % save the included channel names in ch_include

elec_indx = strcmpi(elec_include,'yes');
dataBase(subj).metadata.elec_include = elec_indx;% save a logical vector including all channels and their inclusion yes (1) or no (0) in elec_include for easy computational matters

% check wich channels are included in the stimulation, but excluded by my electrode selection en which are not included in the stimulation but
% included in my electrode selection

stimchannels = unique(dataBase(subj).metadata(run).cc_stimchans); 

j=1;
for i = 1:size(stimchannels,1)
if sum(strcmpi(chan_include,stimchannels{i,1})) == 0
   fprintf('... Channel %s of %s is a bad channel included in stim... \n',stimchannels{i,1},dataBase(subj).sub_label)
   dataBase(subj).metadata.bad_stim{j,1} = stimchannels{i,1};
   j = j+1;
end
end

j = 1;
for i = 1:size(chan_include,1)
if sum(strcmpi(dataBase(subj).metadata(run).cc_stimchans,chan_include{i})) == 0
   fprintf('...Channel %s of %s is not a stimulation channel... \n',chan_include{i},dataBase(subj).sub_label)
   dataBase(subj).metadata.no_stim{j,1} = chan_include{i};
   j = j+1;
end
end
end

clear stimchannels chan_include elec_include elec_indx idx_all idx_bad bad idx_csf idx_screw idx_silicon idx_whitematter j i subjs run
%% mrtrix warning
warning('make sure you have the MRIs in the right directory (see mrtrix code in STReEF02_coreg_roidef_mrtrix and run section 1)')
% run this part after you runned section 1 of 'STReEF02_coreg_roidef_mrtrix'
%% co-register MRI and DWI data with SPM
% co-registration of the MRI acquired at the time of the DWI scan (the 'MRI-DWI') and the DWI. 

for  subj=1:size(dataBase,2)
flags = spm_get_defaults('coreg.estimate');
flags.graphics = ~spm('CmdLine');

sub_label = dataBase(subj).sub_label;

dir = [myDataPath.coreg_ROIpath sprintf('%s/',sub_label)]; % path to the designated co-registration and roi definition folder ('coreg_ROI') 
reference = [dir 'mean_b0_preprocessed.nii'] ; % DWI data as reference image
source = [dir sprintf('MRI_DWI_%s.nii',sub_label)]; % MRI-DWI data as source image

x_dwi = spm_coreg(reference,source,flags); % check the co-registration in the spm window

t_matrix_dwi=spm_matrix(x_dwi(:)'); % extract the transformation matrix
save([dir 'transmatrix_dwi.txt'],'t_matrix_dwi','-ascii') % save the transformation matrix in the designated co-registration and roi definition folder ('coreg_ROI')
dataBase(subj).metadata.transmatrix_dwi = t_matrix_dwi; % save the transformation matrix in the dataBase for inspection and further use
end
clear elec_indx sub_label source reference flags x_dwi

%% SECTION 2:co-register MRI data with SPM, and calculate the transformation matrix 
%% mrtrix warning
warning('make sure you have the transformed MRIs in the right directory (see mrtrix code in STReEF02_coreg_roidef_mrtrix and run section 2)')
% run this part after you runned section 2 of 'STReEF02_coreg_roidef_mrtrix'
%% co-register MRI data with SPM
% co-registration of the MRI acquired before electrode implantation (the 'MRI-CT') and the MRI acquired at the time of the DWI scan (the 'MRI-DWI'). 

for subj = 1:size(dataBase,2)
flags = spm_get_defaults('coreg.estimate');
flags.cost_fun = 'ncc';
flags.graphics = ~spm('CmdLine');

sub_label = dataBase(subj).sub_label;

dir = [myDataPath.coreg_ROIpath sprintf('%s/',sub_label)]; % path to the designated co-registration and roi definition folder ('coreg_ROI') 
reference = [dir sprintf('MRI_DWI_%s_coreg.nii',sub_label)]; % MRI-DWI data as reference image
source = [dir sprintf('%s_ses-1_T1w.nii',sub_label)]; % MRI-CT data as source image

x_mri = spm_coreg(reference,source,flags); %check the co-registration in the spm window

t_matrix_mri=spm_matrix(x_mri(:)');
save([dir 'transmatrix_mri.txt'],'t_matrix_mri','-ascii') % save the transformation matrix in the designated co-registration and roi definition folder ('coreg_ROI')
dataBase(subj).metadata.transmatrix_mri = t_matrix_mri; % save the transformation matrix in the dataBase for inspection and further use
end
clear elec_indx sub_label source reference flags x_mri

%% SECTION 3:transform electrode contact coordinates 
%% mrtrix warning
warning('make sure you have the transformed CTs in the right directory (see mrtrix code in STReEF02_coreg_roidef_mrtrix and run section 3)')
%% transform coordinates
%transform the electrode contact coordinates extracted from the CT scan with the transformation matrix calculated between the 'MRI-CT' and the 'MRI-DWI'. 
% This transformation matrix can be used to bring the coordinates (who are in the MRI-CT space) into the MRI-DWI space. 

for subj=1:size(dataBase,2) 
elec_indx = dataBase(subj).metadata.elec_include;
coordinates = [dataBase(subj).tb_electrodes.x(elec_indx) dataBase(subj).tb_electrodes.y(elec_indx) ...
                dataBase(subj).tb_electrodes.z(elec_indx)]; % the columns are the x,y,z coordinates, the rows the channels

% for some patients the variable tb_electrodes is stored in a different type of array. If so, change.
if isa(coordinates,'cell')
coordinates = str2double(coordinates);
end

% calculate the inverse of the transformation matrix
t_matrix_mri = dataBase(subj).metadata.transmatrix_mri;
rot_trans = [1,0,0,-t_matrix_mri(1,4);0,1,0,-t_matrix_mri(2,4);0,0,1,-t_matrix_mri(3,4);0,0,0,1];
lat_trans = [t_matrix_mri(1:3,1:3)';0,0,0];
lat_trans = [lat_trans,[0;0;0;1]];
t_matrix_inv_mri = lat_trans*rot_trans;

% apply the inversed transformation matrix on the electrode contact coordinates
one = ones(size(coordinates,1),1);
points = [coordinates one]';
coordinates_trans = t_matrix_inv_mri * points; 
coordinates_trans = coordinates_trans(1:3,:)';
dataBase(subj).metadata.coordinates_trans_mri = coordinates_trans;

% make volumes to check if the electrode contact coordinates are transformed correctly
sub_label = dataBase(subj).sub_label;
dir = [myDataPath.coreg_ROIpath sprintf('%s/',sub_label)]; % path to the designated co-registration and roi definition folder ('coreg_ROI') 
CT_BIDS = [dir sprintf('CT_BIDS_%s_coreg.nii',sub_label)]; 
MRI_DWI = [dir sprintf('MRI_DWI_%s_coreg.nii',sub_label)];
DWI = [dir 'mean_b0_preprocessed.nii'];
[output,~,~,outputStruct] = position2reslicedImage(coordinates_trans,MRI_DWI);  % saves the electrode coordinates in the MRI-DWI volume
outputStruct.fname = [dir 'coordinates_trans_mri.nii'];
spm_write_vol(outputStruct,output);
[output,~,~,outputStruct] = position2reslicedImage(coordinates_trans,DWI); % saves the electrode coordinates in the DWI volume
outputStruct.fname = [dir 'coordinates_trans_dwi.nii'];
spm_write_vol(outputStruct,output);
[output,~,~,outputStruct] = position2reslicedImage(coordinates_trans,CT_BIDS); % saves the electrode coordinates in the CT volume
outputStruct.fname = [dir 'coordinates_trans_CT.nii'];
spm_write_vol(outputStruct,output);
end
clear one points lat_trans rot_trans coordinates elec_indx output outputStruct MRI_DWI CT_BIDS DWI t_matrix_dwi t_matrix_inv_mri t_matrix_mri 

%% SECTION 4:define region-of-interests(ROI)/electrode contact areas
%% mrtrix warning
warning('make sure you have the grey-white matter boundary  mask in the right directory (see mrtrix code in STReEF02_coreg_roidef_mrtrix and run section 4)')
%% define electrode contact areas
% make ROI's for every electrode contact, we name those ROI's the electrode contact areas. Make the electrode contact areas by projecting the electrode contact coordinates 
% onto the grey-white matter boundary, define an area by assigning the nearest 64 voxels to each electrode contact area, handle the overlapping voxels, and calculate the volumes of the electrode
% contact areas.

for subj=1:size(dataBase,2)
sub_label = dataBase(subj).sub_label;
dir = [myDataPath.coreg_ROIpath sprintf('%s/',sub_label)]; % path to the designated co-registration and roi definition folder ('coreg_ROI') 
gmwm_mask = read_mrtrix([dir 'gmwmSeed_mask.mif']); % the grey-white matter boundary mask

% extract the coordinates of the grey-white matter boundary mask
logical_mask = logical(gmwm_mask.data); % make the mask binair
check_logical = gmwm_mask;
check_logical.data = logical_mask;
write_mrtrix (check_logical, [dir 'check_logical.mif']) % save the binairy mask for checking the ROI definition steps
[x2, y2, z2] = ind2sub(size(logical_mask),find(logical_mask)); % get the indeces
pos_mask = [x2, y2, z2]; 
vox2real = gmwm_mask.transform; % transform the mask coordinates from the voxel space to the real world space
one = ones(size(pos_mask,1),1);
points = [pos_mask one]';
cor_mask_real = vox2real * points;
cor_mask_real = cor_mask_real(1:3,:)';

coordinates_trans = dataBase(subj).metadata.coordinates_trans_mri; % the electrode contact coordinates in the right MRI-DWI space

% compute the 64 nearest voxels the grey-white matter boundary mask for every electrode contact coordinate and combine to one mask
dir2 = [dir 'coordinate_roi/']; % path to a sub-folder of the designated co-registration and roi definition folder ('coordinate_roi') (made for clarity, not necessary?)
template = [dir 'gmwmSeed_mask.mif']; % path to the grey-white matter boundary mask
roi = zeros(size(gmwm_mask.data));
for cor = 1:size(coordinates_trans,1)
[D2, I2] = pdist2(cor_mask_real,coordinates_trans(cor,:),'euclidean','Smallest',64); % compute the euclidian distance between the electrode contact coordinates and the gray-white matter boundary mask
coordinate_roi= cor_mask_real(I2,:);

% if you want to check this step, use these two lines to save the data
%outname = [dir2 sprintf('coordinate_roi_%d.mif',cor)];
%data = position2reslicedImage_mif(coordinate_roi,template,outname,1);

data = position2reslicedImage_mif_nosave(coordinate_roi,template,1);
roi = roi+data.data;
end
roi_header = data;
roi_header.data = roi;
write_mrtrix (roi_header,[dir2 'coordinate_all.mif']) % save the nearest 64 voxels for every electrode contact area in one volume

% handle the overlapping voxels of the preliminary electrode contact areas
[x4, y4, z4] = ind2sub(size(roi),find(roi > 1)); % % compute the overlap between the preliminary electrode contact areas
pos_overlap = [x4, y4, z4]; 
one = ones(size(pos_overlap,1),1);
points = [pos_overlap one]';
cor_overlap_real = vox2real * points; % transform the overlapping coordinates in the electrode contact areas from the voxel space to the real world space
cor_overlap_real = cor_overlap_real(1:3,:)';
[D4, I4] = pdist2(coordinates_trans,cor_overlap_real,'euclidean','Smallest',1); % compute the euclidian distance between the electrode contact coordinates and the overlapping voxel coordinates in the electrode contact areas
% vector I4 contains the numbers of the closest electorde contact coordinate

% remove per electrode contact area the voxels with overlap and assign each overlap voxel to the closest electrode contact area
volume_roi= NaN(size(coordinates_trans,1),1); 
roi_connectome = zeros(size(gmwm_mask.data));

for cor = 1:size(coordinates_trans,1)
[D2, I2] = pdist2(cor_mask_real,coordinates_trans(cor,:),'euclidean','Smallest',64); % compute the euclidian distance between the electrode contact coordinates and the gray-white matter boundary mask
coordinate_roi = cor_mask_real(I2,:);
[val,pos]=intersect(coordinate_roi,cor_overlap_real,'rows');  % voxels with overlap
coordinate_roi(pos,:)=[]; % remove voxels with overlap
if intersect(I4,cor)
closest = ismember(I4,cor); % see if this electrode contact is in the list of the closest electrode contact area for the each voxel with overlap
coordinate_roi = [coordinate_roi; cor_overlap_real(closest,:)]; % assign the voxel with overlap to the closest electrode contact area
end

data2 = position2reslicedImage_mif_nosave(coordinate_roi,template,cor); % save the coordinates of the electrode contact areas and give those voxels a value, corresponsing to their electrode contact coordinate index
roi_connectome = roi_connectome+data2.data; % store all the coordinates of the final electrode contact areas in one volume

%compute the final volume per electrode contact area (64 mm3 or less)
volume_roi(cor) = size(coordinate_roi,1);

% if you want to check this step, use these two lines to save the data
%outname2 = [dir2 sprintf('coordinate_roi_overlap_%d.mif',cor)];
%data2 = position2reslicedImage_mif(coordinate_roi,template,outname2,cor);
end
% save the output (electrode contact areas + volume) per patient
roi_con_header = data2;
roi_con_header.data = roi_connectome;
write_mrtrix (roi_con_header,[dir2 'coordinate_all_connectome.mif']) % save the final electrode contact areas in one volume
total_volume = sum(volume_roi);
dataBase(subj).metadata.volume_roi = volume_roi;
dataBase(subj).metadata.total_volume = total_volume;
save([dir2 'volume_roi.txt'],'volume_roi','-ascii') % save the volume per electrode contact area
save([dir2 'total_volume.txt'],'total_volume','-ascii') % save the volume of all electrode contact areas combined, per patient
end

clear one points check_logical closest coordinate_roi coordinates_trans cor cor_mask_real cor_overlap_real D2 D4 data data2 dir dir2 gmwm_mask I2 I4 logical_mask pos pos_mask pos_overlap roi roi_con_header roi_connectome roi_header template total_volume val volume_roi vox2real x2 x4 y2 y4 z2 z4
%% SECTION 5: save the important variables of the whole file
%% mrtrix warning
warning('make sure you also run the final section of the twin written in mrtrix code (see STReEF02_coreg_roidef_mrtrix and run section 5)')
%%  save the important variables of the whole file
% put all the database variables in one struct and save it as dwiInfo.mat 
for subj = 1:size(dataBase,2)
    targetFolder = fullfile(myDataPath.DWIMATLABpath,dataBase(subj).sub_label,dataBase(subj).ses_label); % save the computed variables for further use in a designated folder 'dwi_matlab' 

    % Create the folder if it doesn't exist already.
    if ~exist(targetFolder, 'dir')
        mkdir(targetFolder);
    end
    
    start_filename = strfind(dataBase(subj).metadata.dataName,'/');
    stop_filename = strfind(dataBase(subj).metadata.dataName,'_ieeg');
    
    fileName = [dataBase(subj).metadata.dataName(start_filename(end)+1:stop_filename-1),'_dwiInfo.mat']; % name it in a similair way as the CCEP files
    
    % insert variables
    dwi = struct();
    dwi.dataName = dataBase(subj).metadata.dataName;
    dwi.elec_include = dataBase(subj).metadata.elec_include; % logical vector including all channels and their inclusion yes (1) or no (0) for easy computational matters
    dwi.ch = dataBase(subj).metadata.ch; % all channel names
    dwi.cc_stimchans = dataBase(subj).metadata.cc_stimchans; % channel names of stimulated electrode pairs
    dwi.cc_stimsets = dataBase(subj).metadata.cc_stimsets; % channel number of stimulated electrode pairs
    dwi.transmatrix_mri = dataBase(subj).metadata.transmatrix_mri; % transformation matrix between the MRI-CT and the MRI-DWI data
    dwi.transmatrix_dwi = dataBase(subj).metadata.transmatrix_dwi; % transformation matrix between the MRI-DWI and the DWI data
    dwi.coordinates_trans_mri = dataBase(subj).metadata.coordinates_trans_mri; % transformated electrode contact coordinates
    dwi.volume_roi = dataBase(subj).metadata.volume_roi; % volume per electrode contact area

    save(fullfile(targetFolder,fileName), '-struct','dwi');
    fprintf('Saved dwi-struct in %s \n',fullfile(targetFolder,fileName))
end
