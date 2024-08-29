function dataBase = calculate_topology(dataBase)

% INPUT:
% SC_matrix = network modality one
% EC_matrix = network modality two
% elec_info = information about the electrodes (tb_electrodes in BIDS)
% containing channel names and soz information
% elec_include = included electrodes

% OUTPUT:
%

% calculate degree, betweenness centrality, node proximity, and volume of electrode areas (VEA) from the effective and structural networks

for nSubj = 1:size(dataBase,2)

    ch = dataBase(nSubj).ch_select;
    SC_matrix = dataBase(nSubj).SC_matrix;
    EC_matrix = dataBase(nSubj).EC_matrix;

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
    coordinates = [dataBase(nSubj).x_select dataBase(nSubj).y_select ...
        dataBase(nSubj).z_select]; % [channels * x,y,z]

    % compute distance between coordinates
    distances = NaN(size(coordinates,1)); % distance in mm between two electrode contacts
    for  nChan = 1:size(coordinates,1)
        for  mChan = 1:size(coordinates,1)
            D_dis = pdist2(coordinates(nChan,:),coordinates(mChan,:),'euclidean'); % compute the euclidian distance between the electrode contact coordinates
            if ~isempty(D_dis)
                distances(nChan,mChan) = D_dis;
            end
        end
    end

    % compute the node proximity:  the node proximity per node was defined as the median distance between this node and all other nodes.
    node_proximity = NaN(size(ch,1),1);
    for nChan = 1:size(ch,1)
        node_proximity(nChan) = median(nonzeros(distances(:,nChan)),'omitnan');
    end

    % make a table with all the topology measures
    dataBase(nSubj).topology = table(degree_SC,degree_EC,BC_SC,BC_EC,node_proximity);
end

end