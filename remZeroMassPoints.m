function [y,net_edges,y_net_edge_inds] = remZeroMassPoints(y,net_edges,edges,edge_weights,y_net_edge_inds)
%remZeroPointsMass removes points on y with zero projected mass
%   differs from public version in that it takes in FIX_ENDPTS parameter,
%   so that the fixed endpoints don't get removed.

[m,~] = size(y);
y_mass = zeros(m,1); %the mass going through y_i outside the network
net_bool = ismember(edges,fliplr(net_edges),'rows');
non_net_inds = find(~net_bool);

for i=1:m
    y_mass(i) = y_mass(i) + sum(edge_weights(non_net_inds(edges(non_net_inds,2)==i | edges(non_net_inds,1)==i)));
end
removals = y_mass == 0;
removals = sort(removals,'descend'); %%%%%%% EDIT SO THAT TOPOLOGICAL VERTICES DON'T GET REMOVED (function needs extra arguments)
y_new = y;

for i=1:length(removals)
    rem_ind = removals(i);
    y_new(rem_ind,:) = [];
    
    vert_indices = vert_indices - double(vert_indices>=rem_ind);
    
    %edges(edge_inds(rem_ind),:) = [];
    rem_e_inds = y_net_edge_inds{rem_ind};
    rem_edges = net_edges(rem_e_inds,:);
    neigh1 = setdiff(rem_edges(1,:),rem_ind);
    neigh2 = setdiff(rem_edges(2,:),rem_ind);
    %neighs = setdiff(union(net_edges(net_edge_inds{rem_ind},:),[]),rem_ind);
    net_edges = [net_edges;[neigh1,neigh2]];
    net_edges(rem_e_inds,:) = [];
    net_edges = net_edges - double(net_edges>=rem_ind);
    
    y_net_edge_inds{neigh1} = setdiff(y_net_edge_inds{neigh1},rem_e_inds(1));
    y_net_edge_inds{neigh2} = setdiff(y_net_edge_inds{neigh2},rem_e_inds(2));
    for j=1:length(y_net_edge_inds)
        y_net_edge_inds{j} = y_net_edge_inds{j} - double(y_net_edge_inds{j}>rem_e_inds(1))-double(y_net_edge_inds{j}>rem_e_inds(2));
    end
    y_net_edge_inds{neigh1} = [y_net_edge_inds{neigh1},length(net_edges)];
    y_net_edge_inds{neigh2} = [y_net_edge_inds{neigh2},length(net_edges)];
    y_net_edge_inds = {y_net_edge_inds{1:rem_ind-1},y_net_edge_inds{rem_ind+1:end}};
    
end
end

