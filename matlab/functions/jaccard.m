function [JI JI_expected JI_parker JI_alternative] = jaccard(SC_matrix,EC_matrix,dir)
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
% 
% Calculate the inter-modal similarity between two networks

save([dir '/SC_matrix_dwi'],'SC_matrix','-ascii') % save the structural network
save([dir '/EC_matrix_ccep'],'EC_matrix','-ascii') % save the effective network
JI = sum(sum(SC_matrix & EC_matrix))/sum(sum(SC_matrix | EC_matrix)); % calculate the Jaccard Index to be able to compare it with the Jaccard Index calculated with R

% calculate expected Jaccard to check with R
pis= size(EC_matrix(EC_matrix==1),1)/(size(EC_matrix,1)*size(EC_matrix,1)); % calculation how it matches the R calculation, do we agree with this?
pjs = size(SC_matrix(SC_matrix==1),1)/(size(SC_matrix,1)*size(SC_matrix,1));
JI_expected = (pis*pjs)/((pis+pjs)-(pis*pjs)); 

% how I now think we should calculate it:
n_EC = size(EC_matrix,1); % number of nodes
p_EC= (n_EC*(n_EC-1))/2; % potential connections is (n*(n-1))/2)
a_EC = size(EC_matrix(EC_matrix==1),1);

n_SC = size(SC_matrix,1); % number of nodes
p_SC= (n_SC*(n_SC-1))/2; % potential connections is (n*(n-1))/2)
a_SC = size(SC_matrix(SC_matrix==1),1); % actual connections

d_EC = a_EC/p_EC ;  % density is actual connections/potential connections
d_SC = a_SC/p_SC ;  % density is actual connections/potential connections

JI_alternative = (d_EC*d_SC)/(d_SC + d_EC -(d_EC*d_SC)); % this results in very high expected JI's
JI_parker = (d_EC*d_SC)/(d_EC+d_SC); % this is how Christopher Parker explained it to me via email (https://www.sciencedirect.com/science/article/pii/S221315821730325X?via%3Dihub)

end