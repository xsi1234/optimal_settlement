function [y,neighbors,edges,edge_costs] = growVerts(x,y,net_edges,edges,edge_weights,edge_costs,alpha,lambda)
%growVerts
% abandoned function whose intention was to check when it is advantageous
% add a leaf (vertex of one degree) to the network

m = length(y(:,1));
n = length(x(:,1));

x_edges = edges(:,1)>=n;
for i=1:m
    xy_i_edges = x_edges & edges(:,2)==i;
    X_i = x(edges(xy_i_edges,1),:);

end

