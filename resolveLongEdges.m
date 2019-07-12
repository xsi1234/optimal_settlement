function [y,net_edges] = resolveLongEdges(y,net_edges,max_leng)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
m = length(y(:,1));
edge_lengs = sqrt(sum((y(net_edges(:,1),:)-y(net_edges(:,2),:)).^2,2));
long_inds = edge_lengs > max_leng;

new_pts = .5*y(net_edges(long_inds,1),:) + .5*y(net_edges(long_inds,2),:);
num_new_pts = length(new_pts(:,1));
new_edges1 = [(m+1:m+num_new_pts)',net_edges(long_inds,1)];
new_edges2 = [(m+1:m+num_new_pts)',net_edges(long_inds,2)];

y = [y;new_pts];
net_edges(long_inds,:) = [];
net_edges = [net_edges;new_edges1;new_edges2];

end

