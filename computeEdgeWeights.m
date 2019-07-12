function [edges,edge_costs,edge_weights,net_edges,y_net_edge_inds,neighbors,y,z,b] = computeEdgeWeights(y,x,alpha,lambda,net_edges,mass,check_top)
% future optimizations:
% save n^2 ops on saving cost matrix for data

%[N,~] = size(cost_mat); % number of points on graph
n = length(mass); %number of data points
%m = N-n; %number of points on network
[m,~] = size(y);
N = n+m;

% Compute (m+n)x(m+n) cost matrix
dist_mat = squareform(pdist([y;x])); %note that points on y are listed first
net_inds1 = sub2ind(size(dist_mat),net_edges(:,1),net_edges(:,2));
net_inds2 = sub2ind(size(dist_mat),net_edges(:,2),net_edges(:,1));
net_inds = [net_inds1;net_inds2];
dist_mat(net_inds) = alpha*dist_mat(net_inds);
%dist_mat(net_edges(:,1),net_edges(:,2)) = alpha*dist_mat(net_edges(:,1),net_edges(:,2));
%dist_mat(net_edges(:,2),net_edges(:,1)) = alpha*dist_mat(net_edges(:,2),net_edges(:,1));
% Compute all pairs shortest paths for points on network
[dists,next] = floyd(m,dist_mat);
% Compute all pairs shortest paths
% dists = zeros(N,N);
% next = zeros(N,N);
% dists(1:m,1:m) = y_dists;
% next(1:m,1:m) = y_next;
% for i=1:m
%     for j=1:n
%         j_dist_to_net = sqrt(sum((repmat(x(j,:),m,1)-y).^2,2));
%         i_intern_dist = sqrt(sum((repmat(y(i,:),m,1)-y).^2,2));
%         tot_dist = j_dist_to_net + i_intern_dist;
%         [min_dist,min_ind] = min(tot_dist);
%         straight_dist = sqrt(sum((y(i,:)-x(j,:)).^2,2));
%         if min_dist<straight_dist
%             dists(i,m+j) = min_dist;
%             dists(m+j,i) = min_dist;
%             next(i,m+j) = min_ind;
%             next(m+j,i) = min_ind;
%         else
%             dists(i,m+j) = straight_dist;
%             dists(m+j,i) = sraight_dist;
%             next(i,m+j) = m+j;
%             next(m+j,i) = i;
%         end
%     end
% end
%
% for i=1:n
%     for j=i+1:n
%         i_dist_to_net = sqrt(sum((repmat(x(i,:),m,1)-y).^2,2));
%         j_dist_to_net = sqrt(sum((repmat(x(j,:),m,1)-y).^2,2));
%         net_dist = i_dist_to_net + j_dist_to_net;
%         [min_dist,min_ind] = min(net_dist);
%         straight_dist = sqrt(sum((x(i,:)-x(j,:)).^2,2));
%         if min_dist<straight_dist
%             dists(m+i,m+j) = min_dist;
%             dists(m+j,m+i) = min_dist;
%             next(m+i,m+j) = min_ind;
%             next(m+j,m+i) = min_ind;
%         else
%             dists(m+i,m+j) = straight_dist;
%             dists(m+j,m+i) = sraight_dist;
%             next(m+i,m+j) = m+j;
%             next(m+j,m+i) = m+i;
%         end
%     end
% end
% Compute edge weights
weights_mat = zeros(N,N);
for i=m+1:N
    for j=i+1:N
        k = i;
        while next(k,j)>0 %increment weights on edges from shortest path from i to j (and j to i)
            weights_mat(k,next(k,j)) = weights_mat(k,next(k,j)) + 2*mass(i-m)*mass(j-m);
            weights_mat(next(k,j),k) = weights_mat(k,next(k,j));
            k = next(k,j);
        end
    end
end

% for i=1:length(net_edges)
%     weights_mat(net_edges(i,1),net_edges(i,2)) = weights_mat(net_edges(i,1),net_edges(i,2))*alpha + lambda;
%     weights_mat(net_edges(i,2),net_edges(i,1)) = weights_mat(net_edges(i,1),net_edges(i,2));
% end
costs_mat = weights_mat;  %should update costs if net_edges are changed
for i=1:length(net_edges(:,1))
    costs_mat(net_edges(i,1),net_edges(i,2)) = costs_mat(net_edges(i,1),net_edges(i,2))*alpha + lambda;
    costs_mat(net_edges(i,2),net_edges(i,1)) = costs_mat(net_edges(i,1),net_edges(i,2));
end

%k = find(triu(weights_mat)>0);
k = find(triu(costs_mat)>0);
edge_weights = weights_mat(k);
edge_costs = costs_mat(k);
[row,col] = ind2sub(size(weights_mat),k); %the entries in col are increasing
edges = [col,row];
edges(col<row,:) = [row(col<row),col(col<row)]; %first index is always greater than second index
net_edge_inds = edge_weights>lambda/(1-alpha);
if check_top
    net_edges = edges(net_edge_inds,:);
    net_edges = fliplr(net_edges);
    net_edges = net_edges(net_edges(:,2)<=m,:);
end
%num_intern_edges = sum(col<=m); %the internal edges are the first listed in edges_list

removals = setdiff(1:m,union(edges(:,1),edges(:,2)));
for i=1:length(removals)
    edges = edges - double(edges>removals(end-i+1));
    net_edges = net_edges - double(net_edges>removals(end-i+1));
end
y(removals,:) = [];
m = length(y(:,1));
y_net_edge_inds = cell(m,1);
neighbors = cell(m,1);
for i=1:length(net_edges(:,1))
    neighbors{net_edges(i,1)} = sort([neighbors{net_edges(i,1)}, net_edges(i,2)]);
    neighbors{net_edges(i,2)} = sort([neighbors{net_edges(i,2)}, net_edges(i,1)]);
    y_net_edge_inds{net_edges(i,1)} = sort([y_net_edge_inds{net_edges(i,1)},i]);
    y_net_edge_inds{net_edges(i,2)} = sort([y_net_edge_inds{net_edges(i,2)},i]);
end

removals2 = setdiff(1:m,edges((edges(:,1)>m),2)); %remove points where data don't 'go'
removals2 = fliplr(removals2); %so it is in descending order (setdiff sorts them)
%y(removals2,:) = [];
%net_edge_rem_inds = [];
for i=1:length(removals2)
    rem_i = removals2(i);
    neighs = neighbors{rem_i};
    
    if length(neighs) < 3  %only remove points of degree < 3
        y(rem_i,:) = [];
        if ~isempty(neighs)
            neighbors{neighs(1)} = setdiff(neighbors{neighs(1)},rem_i);
            y_net_edge_inds{neighs(1)} = setdiff(y_net_edge_inds{neighs(1)},y_net_edge_inds{rem_i});
            if isempty(neighbors{neighs(1)})
                neighbors{neighs(1)} = []; %just for the assert statement below
                y_net_edge_inds{neighs(1)} = [];
            end
            
            if length(neighs) == 2
                coef = 0;
                if ~ismember(neighs(1),neighbors{neighs(2)}) % then edge will be added to network
                    coef = 1;
                end
                eind1 = edges(:,1) == max(rem_i,neighs(1)) & edges(:,2) == min(rem_i,neighs(1));
                eind2 = edges(:,1) == max(rem_i,neighs(2)) & edges(:,2) == min(rem_i,neighs(2));
                new_weight = edge_weights(eind1);
                assert(new_weight == edge_weights(eind2));
                new_edge = edges(:,1) == max(neighs(1),neighs(2)) & edges(:,2) == min(neighs(1),neighs(2));
                if any(new_edge) % if edge already exists
                    edge_weights(new_edge) = edge_weights(new_edge) + new_weight;
                    edge_costs(new_edge) = edge_costs(new_edge) + new_weight*alpha;
                else % if edge doesn't exist, add it
                    edges = [edges; sort([neighs(1),neighs(2)],'descend')];
                    edge_weights = [edge_weights;new_weight];
                    edge_costs = [edge_costs;new_weight*alpha+lambda*coef];
                end
                
                net_edges = [net_edges;[neighs(1),neighs(2)]];
                neighbors{neighs(2)} = union(setdiff(neighbors{neighs(2)},rem_i),neighs(1));
                neighbors{neighs(1)} = union(neighbors{neighs(1)},neighs(2));
                y_net_edge_inds{neighs(1)} = union(y_net_edge_inds{neighs(1)},length(net_edges));
                y_net_edge_inds{neighs(2)} = union(setdiff(y_net_edge_inds{neighs(2)},y_net_edge_inds{rem_i}),length(net_edges));
            end
            
            rem_edges = edges(:,1) == rem_i | edges(:,2) == rem_i;
            edges(rem_edges,:) = [];
            edge_weights(rem_edges,:) = [];
            edge_costs(rem_edges,:) = [];
            
            rem_net_edges = find(net_edges(:,1) == rem_i | net_edges(:,2) == rem_i);
            net_edges(rem_net_edges,:) = [];
            rev_sorts = sort(rem_net_edges,'descend');
            for j=1:length(y_net_edge_inds)
                    for k=1:length(rev_sorts)
                        y_net_edge_inds{j} = y_net_edge_inds{j} - double(y_net_edge_inds{j}>rev_sorts(k));
                    end
            end
            
%             if ~ismember(y_net_edge_inds{rem_i},net_edge_rem_inds) 
%                 net_edge_rem_inds = [net_edge_rem_inds, y_net_edge_inds{rem_i}];
%                 net_edges(y_net_edge_inds{rem_i},:) = [];
%                 rev_sorts = sort(y_net_edge_inds{rem_i},'descend');
%                 for j=1:length(y_net_edge_inds)
%                     for k=1:length(rev_sorts)
%                         y_net_edge_inds{j} = y_net_edge_inds{j} - double(y_net_edge_inds{j}>rev_sorts(k));
%                     end
%                 end
%             end
            
            
        end
        net_edges = net_edges - double(net_edges>rem_i);
        edges = edges - double(edges>rem_i);
        y_net_edge_inds = {y_net_edge_inds{1:rem_i-1},y_net_edge_inds{rem_i+1:end}}';
        neighbors = {neighbors{1:rem_i-1},neighbors{rem_i+1:end}}';
        for j=1:length(neighbors)
            neighbors{j} = neighbors{j} - double(neighbors{j}>rem_i);
        end
    end
end
if ~isempty(removals2)
    y_net_edge_inds1 = cell(size(neighbors));
    neighbors1 = cell(size(neighbors));
    for i=1:length(net_edges(:,1))
        neighbors1{net_edges(i,1)} = sort([neighbors1{net_edges(i,1)}, net_edges(i,2)]);
        neighbors1{net_edges(i,2)} = sort([neighbors1{net_edges(i,2)}, net_edges(i,1)]);
        y_net_edge_inds1{net_edges(i,1)} = sort([y_net_edge_inds1{net_edges(i,1)},i]);
        y_net_edge_inds1{net_edges(i,2)} = sort([y_net_edge_inds1{net_edges(i,2)},i]);
    end
    assert(isequal(neighbors,neighbors1));
    assert(isequal(y_net_edge_inds,y_net_edge_inds1));
end
% removals = cellfun('isempty', neighbors);
[~,sort_I] = sort(edges(:,1));
edges = edges(sort_I,:);
edge_weights = edge_weights(sort_I);
edge_costs = edge_costs(sort_I);
z = [];
b = [];

end

function [ dists,next ] = floyd(m,weight_mat )
%FLOYD Computes all pairs shortest paths on a graph
%   Edge weights used are given as symmetric matrix weight_mat
% returns: next - a matrix with the (i,j) element specifying the index of
% the first node on the shortest path from i to j.
% dists - a matrix with the (i,j) element specifying the cost of the
% shortest path from i to j. Complexity is O(v^3) where v is the number of
% vertices.

N = length(weight_mat(1,:)); %N = m+n
dists = weight_mat;
next = repmat(1:N,N,1) - diag(1:N);

%for k=1:N
for k=1:m
    for i=1:N
        for j=1:N
            if dists(i,k) + dists(k,j) < dists(i,j)
                dists(i,j) = dists(i,k) + dists(k,j);
                next(i,j) = next(i,k);
            end
        end
    end
end
% adj_indices = dists(m+1:N,m+1:N)>weight_mat(m+1:N,m+1:N);
% [I,J] = ind2sub(size(adj_indices),find(adj_indices));
% dists(adj_indices) = weight_mat(m+I,m+J);
% next(adj_indices) = m+J;
end

function p = path(next,i,j)
p = zeros(1,length(next(:,1)));
k = 1;
p(k) = i;
while next(i,j)>0
    i = next(i,j);
    p(k) = i;
    k = k+1;
end
p(k:end) = [];
end

