function [ y_new,z_new,b_new] = updateyzb(y,z,b,x,rho,edges,edge_costs,first_iter,FIX_ENDPTS)
%UPDATEYZB
% This function performs one step of the ADMM/split bregman algorithm, updating
% each of y,z,b once. 

% INPUTS:
% y - current points for network
% z - variable for ADMM
% b - variable for ADMM
% x - the data points
% rho - coefficient of the penalty for the discrepancy between phi(y) and z
% cut_indices - the last y indices of each component
% edges - matrix with two columns, each representing indices of oriented
%   edges. The first entry of row i represents the starting vertex of edge 
%   i, and we assume edges(i,1)>=edges(i,2) for all i.
% num_intern_edges - the number of internal edges (edges between points on the
%   network). Edge i for i> intern_edges corresponds to an edge between the
%   data and the network.
% edge_costs - vector with entry i corresponding to weight of edge i.
% first_iter - a boolean indicating whether this is the first iteration of
%   ADMM after computing edge weights.
% FIX_ENDPTS - currently needs to be set to 0. Code can be modified below
%   so that if set to 1, ADMM updates keep endpoints of the network fixed.
%
% OUTPUTS:
% y_new - new points for the network
% z_new - new z (dual) variable for ADMM
% b_new - new b variable for ADMM

[m,d] = size(y);
num_intern_edges = sum(edges(:,1)<=m); % edges(:,1) <= edges(:,2) elementwise ALSO ASSUMES edges(:,1) is sorted

assert(sum(edges(:,2)>=edges(:,1))==0);
rel_edge_ind = edges(:,2)<=m;
edges = edges(rel_edge_ind,:); %disregard edges between data points
edge_costs = edge_costs(rel_edge_ind,:);
num_e = length(edges(:,1));
%intern_edges = length(edges(edges(:,1)<=m,:))
z0_1 = y(edges(1:num_intern_edges,2),:) - y(edges(1:num_intern_edges,1),:);
z0_2 = y(edges(num_intern_edges+1:num_e,2),:); %assumes y indices are second coordinate in edges(i)%
c = [zeros(num_intern_edges,d);-x(edges(num_intern_edges+1:num_e,1)-m,:)];
z0 = [z0_1;z0_2] + c;


%Step 1: update z's and b's
if first_iter && ~isempty(b)
    %z0(1:num_intern_edges,:) = z(1:num_intern_edges,:);
    b0 = zeros(size(z0));
    b0(1:num_intern_edges,:) = b(1:num_intern_edges,:);
    b = b0;
else 
    if isempty(b)
        b = zeros(size(z0));
    end
end

prez = z0+b;
normsprez = sqrt(sum(prez.^2,2)); % computes the euclidean norm of each row vector
normsprez(normsprez==0)=1;
normalprez = bsxfun(@rdivide,prez, normsprez);
z = bsxfun(@times,normalprez, max(normsprez-edge_costs/rho,0));

b = b + z0 - z;

z_new = z; b_new = b;

y_new = updateY(y,z_new,b_new,rho,c,edges,num_intern_edges,FIX_ENDPTS);

end




function y_new = updateY(y,z,b,rho,c,edges,intern_edges,FIX_ENDPTS)
%Helper function used for updating y on a single component.

[m,d] = size(y);

A = zeros(m,m);
B = zeros(m,d);
for i=1:length(edges(:,1))
    A(edges(i,2),edges(i,2)) = A(edges(i,2),edges(i,2)) + 1;
    %B(edges(i,2),:) = B(edges(i,2),:) + z(edges(i,2),:) - b(edges(i,2),:) - c(edges(i,2),:);
    B(edges(i,2),:) = B(edges(i,2),:) + z(i,:) - b(i,:) - c(i,:);
    if i<=intern_edges
        A(edges(i,1),edges(i,1)) = A(edges(i,1),edges(i,1)) + 1;
        A(edges(i,1),edges(i,2)) = A(edges(i,1),edges(i,2)) - 1;
        A(edges(i,2),edges(i,1)) = A(edges(i,2),edges(i,1)) - 1;
        %B(edges(i,1),:) = B(edges(i,1),:) - z(edges(i,1),:) + b(edges(i,1),:) + c(edges(i,2),:);
        B(edges(i,1),:) = B(edges(i,1),:) - z(i,:) + b(i,:) + c(i,:);
    end
end


if ~FIX_ENDPTS
    y_new = linsolve(A,B);
else % If endpoints are fixed: %%%% DOES NOT WORK: THIS PART NEEDS TO BE CORRECTED 
    if m == 2
        y_new = y;
    else
        y_new(2:m-1,:) = tridiag((-rho)*ones(m-3,1),D(2:m-1),(-rho)*ones(m-3,1),...
            RHS(2:m-1,:) + rho*[y(1,:);zeros(m-4,d);y(m,:)]);
        y_new([1,m],:) = y([1,m],:);
    end
end

end



