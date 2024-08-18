function topology = calculate_topology(SC_matrix,EC_matrix,elec_info,elec_include)

% INPUT:
% SC_matrix = network modality one
% EC_matrix = network modality two
% elec_info = information about the electrodes (tb_electrodes in BIDS)
% containing channel names and soz information
% elec_include = included electrodes

% OUTPUT:
% 

% calculate degree, betweenness centrality, node proximity, and volume of electrode areas (VEA) from the effective and structural networks

ch_all = elec_info.name; % channel names
ch = ch_all(elec_include);

% graph structure
G_SC = graph(SC_matrix,ch); 
G_EC = graph(EC_matrix,ch);

% degree
degree_SC = degree(G_SC);
degree_EC = degree(G_EC);

% betweennes centrality
BC_SC = centrality(G_SC,'betweenness');
BC_EC = centrality(G_EC,'betweenness');

% node proximity
%extract coordinates
coordinates = [elec_info.x(elec_include) elec_info.y(elec_include) ...
                elec_info.z(elec_include)]; % the columns are the x,y,z coordinates, the rows the channels

% for some patients the variable elec_info is stored in a different type of array. If so, change.
if isa(coordinates,'cell')
coordinates = str2double(coordinates);
end

% compute distance between coordinates
distances = NaN(size(coordinates,1)); % distance in mm between two electrode contacts
for  cor = 1:size(coordinates,1)
for  cor2 = 1:size(coordinates,1)
D_dis = pdist2(coordinates(cor,:),coordinates(cor2,:),'euclidean'); % compute the euclidian distance between the electrode contact coordinates
if ~isempty(D_dis)
    distances(cor,cor2) = D_dis;
else 
    fprintf('empty coordinate in %d or %d',cor,cor2)
end
end
end

% compute the node proximity:  the node proximity per node was defined as the median distance between this node and all other nodes. 
node_proximity = NaN(size(ch,1),1);
for c = 1:size(ch,1)
node_proximity(c) = median(nonzeros(distances(:,c)));
end

% make a table with all the topology measures
topology = table(degree_SC,degree_EC,BC_SC,BC_EC,node_proximity);
end