function [ y,net_edges ] = rndInitialize(x,k)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n = length(x(:,1));
k = min(k,floor(n/2));
%x_k = datasample(x,k,'Replace',false)
ind = union(unidrnd(n,k,1),[]);
rest_ind = setdiff(1:n,ind);
I = dsearchn(x(rest_ind,:),x(ind,:));
repititions = ismember(ind,I);
%new_ind = ind(~repititions);
y = [x(ind(~repititions),:);x(rest_ind(I(~repititions)),:)];
m = length(y(:,1));
net_edges = [(1:m/2)',m/2+(1:m/2)'];


end

