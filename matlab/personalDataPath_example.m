
%% THIS IS AN EXAMPLE FILE
% copy this one and fill in your own

function localDataPath = personalDataPath_example(varargin)

% function that contains local data path, is ignored in .gitignore

if ~isempty(varargin{1})
    if isstruct(varargin{1})
        if sum(contains(fieldnames(varargin{1}),'folderinput'))              
                localDataPath.dataPath = '/blabla/derivatives/STReEF/DATA_for_publication/'; % main data to perform calculations and make figures of manuscript
        end
    end
end
end


end

