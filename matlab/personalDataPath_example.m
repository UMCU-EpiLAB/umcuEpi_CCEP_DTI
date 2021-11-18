
%% THIS IS AN EXAMPLE FILE
% copy this one and fill in your own

function localDataPath = personalDataPath_example(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin{1})
    if isstruct(varargin{1})
        if sum(contains(fieldnames(varargin{1}),'folderinput'))
            if contains(varargin{1}.folderinput,'chronic_ECoG')
                % RESPect database
                localDataPath.CCEPpath = 'blabla/derivatives/CCEP/';
                localDataPath.CCEPpath2 = 'blabla/derivatives/CCEP/'; % second observer
                localDataPath.dataPath = '/blabla/chronic_ECoG/';
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

