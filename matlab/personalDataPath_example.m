
%% THIS IS AN EXAMPLE FILE
% copy this one and fill in your own

function localDataPath = personalDataPath_example(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin{1})
    if isstruct(varargin{1})
        if sum(contains(fieldnames(varargin{1}),'folderinput'))
            if contains(varargin{1}.folderinput,'chronic_ECoG')
                % CCEP
                localDataPath.CCEPpath = 'blabla/derivatives/CCEP/'; % path to the visual scored CCEP data
                localDataPath.CCEPpath2 = 'blabla/derivatives/CCEP/'; % path to the visual scored CCEP data for the second observer (if used)
                localDataPath.dataPath = 'blabla/chronic_ECoG/'; % RESPect database for the ECoGs
                % DTI
                localDataPath.DWIpath = 'blabla/derivatives/preprocess_DWI/' ; % path to the processed DWI's folder
                localDataPath.coreg_ROIpath = 'blabla/derivatives/coreg_ROI/'; % path to the designated co-registration and roi definition folder ('coreg_ROI') 
                localDataPath.DWIMATLABpath = 'blabla/derivatives/dwi_matlab/'; % path to designated folder 'dwi_matlab' to save computed variables in matlab for further use in R/mrtrix 
                localDataPath.FTpath = '/blabla/derivatives/tract2connectome_all/'; % path to the designated fiber tractography folder ('tracts2connectome_all')
            end
        end
    end
end

% % % set paths 
fieldtrip_folder  = '/blabla/fieldtrip/';
% % copy the private folder in fieldtrip to somewhere else
fieldtrip_private = '/blabla/fieldtrip_private/';
addpath(fieldtrip_folder)
addpath(fieldtrip_private)
ft_defaults

end

