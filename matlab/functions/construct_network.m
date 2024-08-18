% construct a symmetric bi-directional effective connectivity network

function dataBase = construct_network(dataBase)

%% construct an effective network matrix
% Connections were drawn from both electrode contacts in a stimulus pair to the electrode contacts in which an ER was detected.  

for subj=1:size(dataBase,2)
    checked = dataBase(subj).ccep.ccep.checked; % matrix with zeros and ones for detected (and visual checked) ERs


    ch = dataBase(subj).ccep.ccep.ch; % all channel names
    elec_include = dataBase(subj).ccep.elec_include; % logical vector including all channels and their inclusion yes (1) or no (0) for easy computational matters
    cc_stimchans = dataBase(subj).ccep.ccep.cc_stimchans; % channel names of stimulated electrode pairs
    ch_include = ch(elec_include); % channel names of the included channels

    EC_matrix_ccep = zeros(size(find(elec_include),1)); % effective matrix with zeros and ones with on the rows and columns the channels

    k=0;
    % run to every stimulated pair to place the detected ERs in the right spot in the matrix with on the rows and columns each included channel one time
    for i = 1:size(find(elec_include),1)
    chan = ch_include(i); % channel name
    loc = strcmpi(chan,cc_stimchans); % stim pairs were this channel was included in
    [col,~] = find(loc); 
    if size(col,1)>1 % if channel was included in two stimulation pairs, combine these results
        col_data = sum(checked(elec_include,col),2); % defines a connection only if there is a ER evoked in after stimulating both stimulation pairs NaN+NaN=NaN, Inf+Inf = Inf Nan+Inf = NaN 1+Inf = NaN 1+NaN = NaN 1+1 = 2, see next line for fix
        [col2,~] = find(checked(elec_include,col)==1); % we want a connection if there is a ER evoked after stimulating one of the stimulation pairs 1+Inf = NaN 1+NaN = NaN 1+1 = 2 assign 1
        col_data(col2,1)=1;
    else
        col_data = checked(elec_include,col); % extract the responses from the included channels after stimulating the selected channel
    end
    
    if ~isempty(col_data)
    EC_matrix_ccep(:,i) = col_data; % insert the response in the effective matrix
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
    end
    end

    EC_matrix_ccep(isnan(EC_matrix_ccep))=0; %make all the NaN zero's
    EC_matrix_ccep(EC_matrix_ccep==Inf)=0; %make al the Inf zero's
    dataBase(subj).ccep.EC_matrix = EC_matrix_ccep; % store the effective matrices in the dataBase for later use
    clear('no_stim_col')

end

disp('EC_matrix made')

%% make the effective network matrix symmetric
% Effective networks were made symmetrical by considering all ERs as bi-directional, to be able to compare to the non-directional structural networks.
for subj=1:size(dataBase,2)
    EC_matrix = dataBase(subj).ccep.EC_matrix;

    EC_matrix_sym = NaN(size(EC_matrix,1));

    for i = 1:size(EC_matrix,1)
    connections = EC_matrix(:,i)|EC_matrix(i,:)'; % consider an connection between channel i and channel j if there is a ER evoked in either direction (so after stimulating channel i or channel j) 
    EC_matrix_sym(:,i) = connections;
    EC_matrix_sym(i,:) = connections;
    
    end
    dataBase(subj).ccep.EC_matrix = EC_matrix_sym;    
end
disp('symmetric matrix made')



end