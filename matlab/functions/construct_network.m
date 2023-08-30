function dataBase = construct_network(dataBase,myDataPath)
% INPUT:
% SC_matrix = network modality one
% EC_matrix = network modality two
% dir = directory to save derivatives

% OUTPUT:
% JI: the jaccard index (JI) which describes the inter-modal similarity between
% two networks
% JI_expected: the expected jaccard index (JI when the connections are
% placed at random positions in the network) as calculated with the R
% function
% JI_parker: the expected jaccard index as explained by Christhoper Parker
% JI_alternative: the expected jaccard index as suggested by me

% construct the effective and structural networks

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
cor_size = size(dataBase(subj).ccep.EC_matrix,2); % is not the same for 5 patients due to one high number in the coordinate_all_connectome.mif file which expanded the matrix with a lot of zeros. I need to figure out why this high number is in the file and fix it in the STReEF02 code.
SC_matrix_dwi = SC_matrix_dwi(1:cor_size,1:cor_size); % for now, workaway by making both the structural and effective networks the same size. Check first if the size difference is really due to the problem described above.
end

dataBase(subj).dwi.SC_matrix_values = SC_matrix_dwi; 

SC_matrix_dwi(SC_matrix_dwi<treshold)=0; % a structural connection was only formed when the streamline density exceeded a threshold of 0.1.
SC_matrix_dwi(SC_matrix_dwi>=treshold)=1; 

dataBase(subj).dwi.SC_matrix = SC_matrix_dwi; 
end
disp('SC_matrix made')

end