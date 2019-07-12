function [ net_edges ] = initializeMST(points)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dist_mat = squareform(pdist(points));
G = graph(dist_mat);
T = minspantree(G);
net_wedges = table2array(T.Edges);
net_edges = net_wedges(:, 1:2);
end

