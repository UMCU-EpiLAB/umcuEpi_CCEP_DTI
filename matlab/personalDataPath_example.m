
%% THIS IS AN EXAMPLE FILE
% copy this one and fill in your own

function localDataPath = personalDataPath_example(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin{1})
    if isstruct(varargin{1})
        if sum(contains(fieldnames(varargin{1}),'folderinput'))  
            if contains(varargin{1}.folderinput,'shareData_STReEF')
            localDataPath.input = '/blabla/STReEF/shareData_STReEF/'; % main effective BIDS data
            localDataPath.input_dev = '/blabla/shareData_STReEF/derivatives/'; % derivatives of the sub scripts, used to run each script on is own  

            % set paths output data
            localDataPath.output = '/blabla/shareData_STReEF/derivatives/'; % save derivatives of the sub scripts    
          
           end
        end   
    end
end

%  set paths fieldtrip 
fieldtrip_folder  = '/blabla/fieldtrip/';
%  copy the private folder in fieldtrip to somewhere else
fieldtrip_private = '/blabla/fieldtrip_private/';
addpath(fieldtrip_folder)
addpath(fieldtrip_private)
ft_defaults


end

