function [JI, JI_expected] = jaccard(SC_matrix,EC_matrix)
% INPUT:
% SC_matrix = network modality one
% EC_matrix = network modality two

% OUTPUT:
% JI: the jaccard index (JI) which describes the inter-modal similarity between
% two networks
% JI_expected: the expected jaccard index (JI when the connections are
% placed at random positions in the network) as calculated with the R
% function

% Calculate the inter-modal similarity between two networks

 % calculate the Jaccard Index 
JI = sum(sum(SC_matrix & EC_matrix))/sum(sum(SC_matrix | EC_matrix));

% calculate the expected Jaccard
pis = size(EC_matrix(EC_matrix==1),1)/(size(EC_matrix,1)*size(EC_matrix,1)); 
pjs = size(SC_matrix(SC_matrix==1),1)/(size(SC_matrix,1)*size(SC_matrix,1));
JI_expected = (pis*pjs)/((pis+pjs)-(pis*pjs)); 

end