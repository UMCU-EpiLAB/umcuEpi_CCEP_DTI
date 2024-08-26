
function localDataPath = STReEF_setLocalDataPath(varargin)

% function LocalDataPath = STReEF_setLocalDataPath(varargin)
% Return the path to the root CCEP  directory and add paths in this repo
%
% input:
%   STReEF_personalDataPath: optional, set to 1 if adding STReEF_personalDataPath
%
% when adding STReEF_personalDataPath, the following function should be in the
% root of this repo:
%
% function localDataPath = STReEF_personalDataPath()
%     'localDataPath = [/my/path/to/data];
%
% this function is ignored in .gitignore
%
% dhermes, 2020, Multimodal Neuroimaging Lab
% dvanblooijs, 2020, UMCU_EpiLAB

if isempty(varargin)

    rootPath = which('STReEF_setLocalDataPath');
    ccepRepoPath = fileparts(rootPath);
    
    % add path to functions
    addpath(genpath(ccepRepoPath));
    
    % add localDataPath default
    localDataPath = fullfile(ccepRepoPath,'data');

elseif ~isempty(varargin)
    % add path to data
    if varargin{1}==1 && exist('STReEF_personalDataPath','file')

        localDataPath = STReEF_personalDataPath();

    elseif varargin{1}==1 && ~exist('STReEF_personalDataPath','file')

        sprintf(['add STReEF_personalDataPath function to add your localDataPath:\n'...
            '\n'...
            'function localDataPath = STReEF_personalDataPath()\n'...
            'localDataPath.input = [/my/path/to/data];\n'...
            'localDataPath.input_dev = [/my/path/to/derivatives];\n'...
            'localDataPath.output = [/my/path/to/output];\n'...
            '\n'...
            'this function is ignored in .gitignore'])
        return
    end

    % add path to functions
    rootPath = which('STReEF_setLocalDataPath');
    ccepRepoPath = fileparts(rootPath);
    addpath(genpath(ccepRepoPath));
    
end

return

