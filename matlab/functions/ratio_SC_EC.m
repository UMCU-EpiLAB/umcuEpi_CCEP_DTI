function ratio = ratio_SC_EC(SC_matrix,EC_matrix)
% INPUT:
% SC_matrix = network modality one
% EC_matrix = network modality two

% OUTPUT:
% Ratio: The ratio between structural and effective connections in the set of symmetric difference connections

% Calculate the inter-modal similarity between two networks

intersecting = SC_matrix - EC_matrix; %set of symmetric difference connections (1 and -1)

ratio = log10(sum(sum(intersecting == 1))/sum(sum(intersecting == -1))); %structural divided by effective
end
