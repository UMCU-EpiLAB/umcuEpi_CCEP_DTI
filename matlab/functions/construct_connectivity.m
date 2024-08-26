% construct a symmetric bi-directional effective connectivity

function dataBase = construct_connectivity(dataBase)

%% construct an effective connectivity matrix
% Connections were drawn from both electrode contacts in a stimulus pair to the electrode contacts in which a CCEP was detected.

for nSubj = 1:size(dataBase,2)

    checked = dataBase(nSubj).ccep.checked; % matrix with zeros and ones for detected (and visual checked) CCEPs
    cc_stimsets = dataBase(nSubj).ccep.cc_stimsets;
    ch = dataBase(nSubj).ccep.ch;
    
    EC_matrix = zeros(size(ch,1),size(ch,1)); % effective connectivity matrix 

    for nStimp = 1:size(cc_stimsets,1)
        EC_matrix(:,cc_stimsets(nStimp,1)) = EC_matrix(:,cc_stimsets(nStimp,1)) + checked(:,nStimp);
        EC_matrix(:,cc_stimsets(nStimp,2)) = EC_matrix(:,cc_stimsets(nStimp,2)) + checked(:,nStimp);
    end

    % make the matrix symmetrical
    EC_matrix = EC_matrix + EC_matrix.';
    EC_matrix(EC_matrix>0) = 1;

    dataBase(nSubj).ccep.EC_matrix = EC_matrix;
end

end