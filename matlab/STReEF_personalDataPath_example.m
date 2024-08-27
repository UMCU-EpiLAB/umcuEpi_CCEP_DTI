
%% THIS IS AN EXAMPLE FILE
% copy this one and fill in your own

function localDataPath = STReEF_personalDataPath_example()

% function that contains local data path, is ignored in .gitignore

localDataPath.input = '/blabla/STReEF/shareData_STReEF/'; % main effective BIDS data
localDataPath.input_dev = '/blabla/shareData_STReEF/derivatives/'; % derivatives of the sub scripts, used to run each script on is own

% set paths output data
localDataPath.output = '/blabla/shareData_STReEF/derivatives/'; % save derivatives of the sub scripts


%  set paths fieldtrip 
fieldtrip_folder  = '/blabla/fieldtrip/';
%  copy the private folder in fieldtrip to somewhere else
fieldtrip_private = '/blabla/fieldtrip_private/';
addpath(fieldtrip_folder)
addpath(fieldtrip_private)
ft_defaults


end

