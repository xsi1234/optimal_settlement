function [edges,edge_costs,edge_weights,net_edges,y_net_edge_inds,neighbors,y,z,b,next,top_change,dists] = computeEdgeWeights(y,x,alpha,lambda,net_edges,mass,...
    check_top,e_heur,delta)
%COMPUTEEDGEWEIGHTS computes edge weights/masses given network
%   This function computes geodesics corresponding to the weighted graph
%   XuY (where the weights correspond to distances, which are multiplied by
%   alpha if the edge lies in the network). 
%OUTPUTS:
% edges - a list of edges on the graph (XuY), whose vertices consist of all
%   data points 'x', and all network points 'y'. This is not the list of
%   the complete edge set, just those edges with positive edge_cost
% edge_costs - a list of costs corresponding to 'edges' (of the graph
%   XuY). These are the energy contributions 'M_{a,b}' defined in the
%   paper: the total energy is \sum |a-b| M_{a,b}.
% edge_weights - a list of weights corresponding to 'edges' (of the graph
%   XuY). This vector contains the total mass traveling through each edge,
%   'm_{a,b}' as defined in the paper.
% net_edges - a list of edges for the initial network. Each edge contains
%   the indices of the points in y0 that form the edge. It is given as a
%   n_e x 2 vector, where n_e is the number of edges
% y_net_edge_inds - a cell that indicates the network edge indices 
%   corresponding to each vertex on y. y_net_edge_inds{i} is a vector of 
%   indices corresponding to the edge number (row) where i appears in 
%   'net_edges'.
% neighbors - a cell that indicates the neighbors of each vertex on the
%   network y. neighbors{i} is a vector of indices which share an edge with
%   vertex i.
% y - locations of points on final network, given as mxd vector
%   (m-number of points, d-dimension). It may change in this function since
%   unnecessary points will be removed.
% z - dual variable of 'y' for ADMM, currently returned as empty.
% b - variable for ADMM, currently returned as empty.
% next - an (m+n)x(m+n) matrix of indices indicating the shortest paths
%   between all pairs of vertices on the graph XuY. the entry at next(i,j)
%   gives the point that should be visited next if going from point i to
%   point j. 
% top_change - value 1 if the topology if the network has, 0 otherwise
% dists - the geodesic (shortest path) distances between all m+n vertices 
%   on the complete graph XuY, with weights equal to length of the edge,
%   times alpha if the edge belongs to the network.
%INPUTS:
% y - locations of points on final network, given as mxd vector
%   (m-number of points, d-dimension)
% net_edges - a list of edges for the initial network. Each edge contains
%   the indices of the points in y0 that form the edge. It is given as a
%   n_e x 2 vector, where n_e is the number of edges
% x - the n data points, given as nxd vector
% mass - the data masses, given as 1xn vector
% lambda - parameter in ONST functional (value > 0)
% alpha - cost of travelling along network, value between 0 and 1.
% check_top - a value equal to 1 if checking to change the topology, 0 
%   otherwise.
% e_heur - a value equal to 1 if checking to change the topology via the
%   heuristic for adding edges (based on ratio of extrinsic to intrinsic 
%   distance), 0 otherwise. To have an effect, check_top should be 1. If
%   check_top = 1 and e_heur = 0, then topology may be changed by simple
%   criteria on how m_{a,b} relates to lambda/(1-alpha).
% delta - the longest distance of potential edges to be added using
%   heurisitc.

%[N,~] = size(cost_mat); % number of points on graph
n = length(mass); %number of data points
%m = N-n; %number of points on network
[m,~] = size(y);
N = n+m;

% Compute (m+n)x(m+n) cost matrix
eucdist_mat = squareform(pdist([y;x])); %note that points on y are listed first
net_inds1 = sub2ind(size(eucdist_mat),net_edges(:,1),net_edges(:,2));
net_inds2 = sub2ind(size(eucdist_mat),net_edges(:,2),net_edges(:,1));
net_inds = [net_inds1;net_inds2];
dist_mat = eucdist_mat;
dist_mat(net_inds) = alpha*eucdist_mat(net_inds);
% Compute all pairs shortest paths for points on network
[edges,edge_costs,edge_weights,net_edges,dists,next,top_change] = getCostsWeights(y,N,dist_mat,net_edges,mass,alpha,lambda,check_top);
energy = calculateEnergy(y,x,edges,edge_costs);
%num_intern_edges = sum(col<=m); %the internal edges are the first listed in edges_list

%if check_top
if check_top && e_heur
    %Check adding edges with large ratio of intrinsic to extrinsic distance
    dist_ratios = triu(dists(1:m,1:m)./eucdist_mat(1:m,1:m),1);
    dist_ratios(eucdist_mat(1:m,1:m)==0) = 0;
    dist_ratios(eucdist_mat(1:m,1:m)>delta) = 0; %don't check vertices further than delta apart
    net_inds1 = sub2ind(size(dist_ratios),net_edges(:,1),net_edges(:,2));
    net_inds2 = sub2ind(size(dist_ratios),net_edges(:,2),net_edges(:,1));
    net_inds = [net_inds1;net_inds2];
    dist_ratios(net_inds) = 0;
    counts = histcounts([net_edges(:,1);net_edges(:,2)],1:m+1);
%     nonvert_inds = counts==2; %only check between topological vertices
%     dist_ratios(nonvert_inds,:) = 0;
%     dist_ratios(:,nonvert_inds) = 0;
    vert_inds = counts~=2; % only check between non-topological vertices
    dist_ratios(vert_inds,:) = 0;
    dist_ratios(:,vert_inds) = 0;
    dist_ratios = dist_ratios.^2;
    if any(any(dist_ratios>0))
        [i,j] = ind2sub([m,m],randsample(m^2,1,true,dist_ratios(1:m^2)));
        if ~edgeIntersects(y,net_edges,[i,j])
            edge_len = eucdist_mat(i,j); %discretize new edge
            net_len = sum(sqrt(sum((y(net_edges(:,1),:)-y(net_edges(:,2),:)).^2,2)));
            avg_elen = net_len/length(net_edges(:,1));
            n_new_e = ceil(edge_len/avg_elen);
            if n_new_e>1
                new_edges = [i,m+1;(m+1:m+n_new_e-2)',(m+2:m+n_new_e-1)';m+n_new_e-1,j];
            else
                new_edges = [i,j];
            end
            net_edges_new = [net_edges;new_edges]; 
            t = (1/n_new_e:1/n_new_e:1-1/n_new_e)';
            new_points = (1-t).*y(i,:) + t.*y(j,:);
            n_new_p = length(new_points(:,1));
            y_new = [y;new_points]; 
            %eucdist_mat and dist_mat need extending to account for new points
            eucdist_mat_new = zeros(N+n_new_p,N+n_new_p);
            old_inds = [1:m,m+n_new_p+1:N+n_new_p];
            eucdist_mat_new(old_inds,old_inds) = eucdist_mat;
            new_dists = pdist2(new_points,[y;x]);
            eucdist_mat_new(m+1:m+n_new_p,old_inds) = new_dists;
            eucdist_mat_new(old_inds,m+1:m+n_new_p) = new_dists';
            eucdist_mat_new(m+1:m+n_new_p,m+1:m+n_new_p) = pdist2(new_points,new_points);
            net_inds1 = sub2ind(size(eucdist_mat_new),net_edges_new(:,1),net_edges_new(:,2));
            net_inds2 = sub2ind(size(eucdist_mat_new),net_edges_new(:,2),net_edges_new(:,1));
            net_inds = [net_inds1;net_inds2];
            dist_mat_new = eucdist_mat_new;
            dist_mat_new(net_inds) = alpha*eucdist_mat_new(net_inds);
            N_new = N + n_new_p;
            [edges_new,edge_costs_new,edge_weights_new,net_edges_new,dists_new,next_new] = getCostsWeights(y_new,N_new,dist_mat_new,net_edges_new,mass,alpha,lambda,0);
            energy_new = calculateEnergy(y_new,x,edges_new,edge_costs_new);
            fprintf('new energy %.8f', energy_new);
            %display(y([i,j],:))
            if energy_new<energy
                edges=edges_new; edge_costs=edge_costs_new; edge_weights=edge_weights_new; net_edges=net_edges_new; dists=dists_new; next = next_new;
                y = y_new; m = length(y_new(:,1));
                fprintf('\n EDGE ADDED')
            end
        end
    end
end

% pos_edges = edges(edge_weights>0,:);
% removals = setdiff(1:m,union(pos_edges(:,1),pos_edges(:,2))); %remove points that have no mass going through them
% for i=1:length(removals)
%     rem_bool = edges(:,1) == removals(end-i+1) | edges(:,2) == removals(end-i+1);
%     edges(rem_bool,:) = [];
%     edge_costs(rem_bool,:) = [];
%     edge_weights(rem_bool,:) = [];
%     rem_bool = net_edges(:,1) == removals(end-i+1) | net_edges(:,2) == removals(end-i+1);
%     net_edges(rem_bool,:) = [];
%     %make the above lines more efficient?
%     edges = edges - double(edges>removals(end-i+1));
%     net_edges = net_edges - double(net_edges>removals(end-i+1));
% end
removals = setdiff(1:m,union(edges(:,1),edges(:,2))); %remove points that are not used
for i=1:length(removals)
    edges = edges - double(edges>removals(end-i+1));
    net_edges = net_edges - double(net_edges>removals(end-i+1));
end
y(removals,:) = [];
m = length(y(:,1));
y_net_edge_inds = cell(m,1);
neighbors = cell(m,1);
for i=1:length(net_edges(:,1))
    neighbors{net_edges(i,1)} = union([neighbors{net_edges(i,1)}], net_edges(i,2));
    neighbors{net_edges(i,2)} = union([neighbors{net_edges(i,2)}], net_edges(i,1));
    y_net_edge_inds{net_edges(i,1)} = sort([y_net_edge_inds{net_edges(i,1)},i]);
    y_net_edge_inds{net_edges(i,2)} = sort([y_net_edge_inds{net_edges(i,2)},i]);
end

[~,sort_I] = sort(edges(:,1));
edges = edges(sort_I,:);
edge_weights = edge_weights(sort_I);
edge_costs = edge_costs(sort_I);
z = [];
b = [];

end



function [edges,edge_costs,edge_weights,net_edges,dists,next,top_change] = getCostsWeights(y,N,dist_mat,net_edges,mass,alpha,lambda,check_top)

m = length(y(:,1));
[dists,next] = floyd(m,dist_mat); %note the matrix 'next' may be inaccurate after any potential topological changes (further below)
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

costs_mat = weights_mat; 
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
cr = col<row;
edges(cr,:) = [row(cr),col(cr)]; %ensures first index is always greater than second index
top_change = 0;
if check_top %add/remove only one edge with highest/lowest weight
    %[max_weight,max_ind] = max(edge_weights);
    flipinds = net_edges(:,1)<net_edges(:,2);
    net_edges(flipinds,:) = [net_edges(flipinds,2),net_edges(flipinds,1)];
    [net_bool,net_locs] = ismember(edges,net_edges,'rows');  %inefficient?
    
    [max_w,max_ind] = max(edge_weights(~net_bool & edges(:,1)<=m,:));  %only adds edges on y (can change to allow adding edges to data)
    if max_w>lambda/(1-alpha)
        net_inds = find(~net_bool & edges(:,1)<=m);
        new_edge = edges(net_inds(max_ind),:);
        if ~edgeIntersects(y,net_edges,new_edge)
            net_edges = [net_edges;new_edge];
            edge_costs(net_inds(max_ind)) = edge_weights(net_inds(max_ind))*alpha + lambda;
            top_change = 1;
        end
    end
    
    zero_inds = edge_weights(net_bool,:)==0;
    if any(zero_inds) %Remove edges with zero mass
        net_loc_inds = net_locs(net_locs>0);
        net_rem_inds = net_loc_inds(zero_inds);
        net_inds = find(net_bool);
        edge_rem_inds = net_inds(zero_inds);
        net_edges(net_rem_inds,:) = [];
        edges(edge_rem_inds,:) = [];
        edge_costs(edge_rem_inds,:) = [];
        edge_weights(edge_rem_inds,:) = [];
        %net_bool(edge_rem_inds) = [];
        %net_locs(edge_rem_inds) = [];
    else
        %Remove edge with smallest weight below lambda/(1-alpha)
        [min_w,min_ind] = min(edge_weights(net_bool,:));
        if min_w<lambda/(1-alpha)
            rem_ind = net_locs(net_locs>0);
            rem_ind = rem_ind(min_ind);
            %net_rem_inds = [net_rem_inds;net_loc_inds(min_ind)];
            net_edges(rem_ind,:) = [];
            net_inds = find(net_bool);
            edge_costs(net_inds(min_ind)) = edge_weights(net_inds(min_ind));
            top_change = 1;
        end
    end

    net_edges = fliplr(net_edges);
      %Add all edges simultaneously
%     net_edge_inds = edge_weights>lambda/(1-alpha);
%     if sum(net_edge_inds>0)
%         net_edges = edges(net_edge_inds,:);
%         net_edges = fliplr(net_edges);
%         net_edges = net_edges(net_edges(:,2)<=m,:);
%         costs_mat = weights_mat;  %should update costs if net_edges are changed
%         for i=1:length(net_edges(:,1))
%             costs_mat(net_edges(i,1),net_edges(i,2)) = costs_mat(net_edges(i,1),net_edges(i,2))*alpha + lambda;
%             costs_mat(net_edges(i,2),net_edges(i,1)) = costs_mat(net_edges(i,1),net_edges(i,2));
%         end
%     end
end
%edge_costs = costs_mat(k);

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



function intersect_bool = edgeIntersects(y,net_edges,new_edge)

ne = length(net_edges(:,1));
x0 = y(new_edge(1),:); 
v0 = y(new_edge(2),:) - y(new_edge(1),:);
intersect_bool = false;
for i=1:ne
    x1 = y(net_edges(i,1),:);
    v1 = y(net_edges(i,2),:) - y(net_edges(i,1),:);
    st = linsolve([-v0',v1'],x0'-x1');
    if all(st<1 & st>0) %check if lines x0 + s*v0, x1 + t*v1 intersect
        intersect_bool = true;
        break;
    end
end
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

