function [ energy ] = calculateEnergy1( y,x,edges,edge_costs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

v = [y;x];
dists = sqrt(sum((v(edges(:,1),:) - v(edges(:,2),:)).^2,2));
energy = dists' * edge_costs;

end

