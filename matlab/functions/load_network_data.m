% function to read derivatives data into dataBase
% do it in the dateMain format (change in dataBase but rest the same!)
% author: Dorien van Blooijs
% date: June 2019

function dataBase = load_network_data(myDataPath,cfg)

input = myDataPath.input;
input_dev = myDataPath.input_dev;

dataBase = struct([]);

for i=1:size(cfg.sub_label,2)
    if isfield(cfg,'sub_label')
        sub_label = ['sub-' cfg.sub_label{i}];
    else
        error('No sub-label specified');
    end
    if isfield(cfg,'ses_label')
        ses_label = cfg.ses_label{i};
    else
        error('No ses-label specified');
    end
        
        
    % load electrodes
    D = dir(fullfile(input,sub_label,ses_label,'ieeg',...
        [sub_label '_' ses_label ,'_electrodes.tsv']));
    
    elecsName = fullfile(D(1).folder, D(1).name);
    
    % remove electrodes named 'other' (so no grid, depth,strip)
    tb_electrodes = readtable(elecsName,'FileType','text','Delimiter','\t');
    idx_elec_incl = ~strcmp(tb_electrodes.group,'other');
    tb_electrodes = tb_electrodes(idx_elec_incl,:);
    
    % load EC and SC derivatives data
    fileName = [sub_label,'_',ses_label,'_Effective_Connectivity.mat'];
    fileName2 = [sub_label,'_',ses_label,'_Structural_Connectivity.mat']; 
    EC = load(fullfile(input_dev,sub_label,fileName));
    SC = load(fullfile(input_dev,sub_label,fileName2));

    dataBase(i).sub_label = sub_label;
    dataBase(i).ses_label = ses_label;
    dataBase(i).tb_electrodes = tb_electrodes;
    dataBase(i).elec_include = EC.elec_include;
    dataBase(i).VEA = SC.VEA;
    dataBase(i).EC_matrix = EC.EC_matrix;
    dataBase(i).SC_matrix = SC.SC_matrix;
    fprintf('...Subject %s has been run...\n',sub_label)
   
end

disp('All subjects are loaded')

