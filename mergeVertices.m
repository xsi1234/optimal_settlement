function [ y, net_edges ] = mergeVertices(y,net_edges,delta)
%UNTITLED Merge points in y, if their distance is less than delta

dists = squareform(pdist(y));

[I_temp,J_temp] = ind2sub(size(dists),find(dists<delta));
if ~isempty(I_temp)
    I = I_temp(I_temp>J_temp);
    J = J_temp(I_temp>J_temp);
    
    [I,sort_ind] = sort(I,'descend');
    J = J(sort_ind);
    if length(I)>1
        I_repeated = I(1:end-1) == I(2:end);
        I(I_repeated) = [];
        J(I_repeated) = [];
    end
    hold on;
    for i=1:length(I)
        %scatter(y(I(i),1),y(I(i),2),100,'black','fill');
        %edges_indices = net_edges(:,1) == I(i) | net_edges(:,2) == I(i);
        %new_edges = net_edges(edges_indices,:);
        net_edges(net_edges(:,1)==I(i),1) = J(i);
        net_edges(net_edges(:,2)==I(i),2) = J(i);
        net_edges = net_edges - double(net_edges>I(i));
        J = J - double(J>I(i));
        
    end
    net_edges(net_edges(:,1)==net_edges(:,2),:) = []; 
    y(I,:) = [];
end

