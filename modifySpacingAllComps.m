function [y_new,z_new,b_new,net_edges,y_net_edge_inds] = modifySpacingAllComps(y,r,top_struct,num_add,max_m,eo,p,FIX_ENDPTS)
% This function makes sure that the spacing between the points on the curve
% is appropriate. If a particular point on the curve does not have
% neighbors on the curve that are close enough, a point is added. Likewise
% if the neighbors are too close, that point is removed. Re-spacing is done
% on all components simultaneosly, and at every function call, either even
% or odd numbered indices are "modified", as specified by input parameter
% 'eo'. Additional points may be added, as indicated by the input parameter
% 'num_add'.
% FUTURE MODIFICATIONS: 1. at the moment only pairs of points are added at a time
% (see calls to addPoints()). It may be more efficient to do this once for all
% points to be added at the same time. Same for removePoints().
% 2. do we want to be 'extending' the endpoints when adding there? This is
% what is currently being done.
%
% INPUTS: y-current points for curve
% mass-the masses of the data points
% I-the assignment of the data points to y (corresponding to shortest distance)
% x-the data points
% r-the threshold for modifying the spacing of points on y
% cut_indices-the last y indices for each component
% num_add-the net gain in number of points on the curve to be desired
% eo-whether to check the odd or even number indices: 0-even, 1-odd
% p: alpha_i^(p-1) * y_mass(i) should be constant
%
% OUTPUTS: y_new - the new points for the curve
% I_new - the estimated assignments for y_new (only the assignments for
% the two neighbors of added/removed points are re-evaluated)
%
%
% Remark: If p=1, it may be desired for points to be added to components of
% size 1 (singletons). The current function doesn't handle this.
%
% Differs from public version in that it uses remZeroMassPoints function
% that has FIX_ENDPTS parameter

assert(r>=1, 'The value of r should be at least 1.');
assert(eo==0 || eo==1, 'The value of eo should be either 0 or 1.'); %FOR NETWORKS EO (EVEN-ODD) DOESN'T WORK THE SAME WAY: NEED FIXING

%[cut_indices,vert_indices,vert_neighs,neighbors,arcs,net_edges,edges,edge_weights,arc_inds,net_edge_inds] = top_struc{:};
[edges,edge_weights,net_edges,y_net_edge_inds,neighbors] = top_struct{:};
vert_indices = find(cellfun(@length,neighbors) > 2);
end_pts_singletons = find(cellfun(@length,neighbors) < 2);
singletons = cellfun(@length,neighbors(end_pts_singletons))<1;
singletons = end_pts_singletons(singletons);
end_pts = setdiff(end_pts_singletons,singletons);

hold on;
% compute mass at each point in y, and remove points with zero mass.
%[y,z_new,b_new,I,y_mass,cut_indices] = remZeroMassPoints(y,z,b,mass,I,cut_indices,FIX_ENDPTS);
[m,d] = size(y);
%edge_inds = cell(m,1); %entry i gives the edges indices that y_i neighbors
y_mass = zeros(m,1); %stores the amount of mass going through each point on y
% for i=1:length(edges)
%     %edge_inds{edges(i,1)} = [edge_inds{edges(i,1)},i];
%     %edge_inds{edges(i,2)} = [edge_inds{edges(i,2)},i];
%     if edges(i,1) <= m
%         y_mass(edges(i,1)) = y_mass(edges(i,1)) + (edge_weights(i) - lambda)/alpha;
%         y_mass(edges(i,2)) = y_mass(edges(i,2)) + (edge_weights(i) - lambda)/alpha;
%     else
%         if edges(i,2) <= m
%             y_mass(edges(i,2)) = y_mass(edges(i,2)) + (edge_weights(i) - lambda)/alpha;
%         end
%     end
% end

%proj_indices = find(edges(:,1)>m & edges(:,2)<=m);  %Change so that it also includes all path weights off network (not only data)
net_bool = ismember(edges,fliplr(net_edges),'rows');
non_net_inds = find(~net_bool);
%proj_indices = find(edges
for i=1:m
    %y_mass(i) = y_mass(i) + sum(edge_weights(proj_indices(edges(proj_indices,2)==i)));
    y_mass(i) = y_mass(i) + sum(edge_weights(non_net_inds(edges(non_net_inds,2)==i | edges(non_net_inds,1)==i)));
end
%y_mass(edges(proj_indices,2)) = edge_weights(proj_indices);


if p==1
    am = y_mass;
else
    % compute alpha(i): the relative length of the interval around y(i)
    % extending to the midpoints of it's 2 neighbors.
    w = zeros(m+2,d);
    w(1,:) = y(1,:);
    w(2:m+1,:) = y(:,:);
    w(m+2,:) = y(m,:);
    
    v = w(2:m+2,:)-w(1:m+1,:);
    normv = sqrt(sum(v.^2,2));
    normv(cut_indices+1) = 0;
    total_length = sum(normv);
    alpha = .5*(total_length/m)*(normv(2:m+1) + normv(1:m));
    if p>2
        alpha = alpha.^(p-1);
    end
    %compute alpha(i)^(p-1)*y_mass(i): want it to be roughly constant for each i
    am = bsxfun(@times,alpha,y_mass);
end
am_med = mean(am);
%am_med = median(am);

%comp_lengs = [cut_indices;m]-[0;cut_indices];
%comp_start_indices = [1;cut_indices+1];
%nonsingleton_indices = setdiff(1:m,comp_start_indices(comp_lengs == 1));

%start_pts = [1;cut_indices+1]; %the start and end indices for all components of y
%end_pts = [cut_indices;m];     %to be checked in the loop below.
%nonendpts = setdiff(1:m,[start_pts;end_pts]);
%nonverts = setdiff(1:m,union(vert_indices,end_pts_singletons));
nonverts = setdiff(1:m,union(vert_indices,end_pts));

[am_sorted,sort_ind] = sort(am);
removals_bool0 = am_sorted==0 & ismember(sort_ind,nonverts);
%removals0 = sort_ind(removals_bool0);
%k0 = length(removals0);
k0 = sum(removals_bool0);
%removals_bool = (am_sorted(1:floor(m/2))/am_med < 1/r) & mod(sort_ind(1:floor(m/2)),2)==eo;
if FIX_ENDPTS
    removals_bool = (am_sorted(1:floor(m/2))/am_med < 1/r) & mod(sort_ind(1:floor(m/2)),2)==eo ...
        & ismember(sort_ind(1:floor(m/2)),nonverts) & ~ismember(sort_ind(1:floor(m/2)),[1,m]);

else
    removals_bool = (am_sorted(1:floor(m/2))/am_med < 1/r) & mod(sort_ind(1:floor(m/2)),2)==eo ...
        & ismember(sort_ind(1:floor(m/2)),nonverts);
end
%removals = sort(sort_ind([removals_bool; boolean(zeros(ceil(m/2),1))]),'descend');
%removals = sort_ind([removals_bool; boolean(zeros(ceil(m/2),1))]);
removals = sort_ind(removals_bool0 | [removals_bool; boolean(zeros(ceil(m/2),1))]); %removals have either deg 2 or 0

k = length(removals);
if ~isempty(num_add) %set k2: the number of points to add
    k2 = min(m-k,k-k0+num_add);
    k2 = min(k2,floor(m/2));
else
    k2 = k;
end
if ~isempty(max_m)
    k2 = max(0,min(max_m-m,k2));
end

%Add points between topological vertices:
top_verts = union(vert_indices,end_pts);
split_edges_inds = find(ismember(net_edges(:,1),top_verts) & ismember(net_edges(:,2),top_verts));
split_edges = net_edges(split_edges_inds,:);
split_num = length(split_edges(:,1));
y_new = [y;zeros(split_num,d)];
ne = length(net_edges(:,1));
net_edges = [net_edges;zeros(split_num,2)];
neighbors_new = [neighbors;cell(split_num,1)];
y_net_edge_inds = [y_net_edge_inds;cell(split_num,1)];

for i=1:split_num
    y_new(m+1,:) = .5*y(split_edges(i,1),:) + .5*y(split_edges(i,2),:);
    net_edges(split_edges_inds(i),:) = [split_edges(i,1),m+1]; %overwrite old edge
    net_edges(ne+1,:) = [split_edges(i,2),m+1]; %add second edge
    y_net_edge_inds{split_edges(i,2)} = [setdiff(y_net_edge_inds{split_edges(i,2)},split_edges_inds(i)),ne+1];
    y_net_edge_inds{m+1} = [split_edges_inds(i),ne+1];
    neighbors_new{split_edges(i,1)} = [setdiff(neighbors_new{split_edges(i,1)},split_edges(i,2)),m+1];
    neighbors_new{split_edges(i,2)} = [setdiff(neighbors_new{split_edges(i,2)},split_edges(i,1)),m+1];
    neighbors_new{m+1} = [split_edges(i,1),split_edges(i,2)];
    m = m + 1;
    ne = ne + 1;
end
k2 = max(0,k2 - split_num);

%split = sort_ind(m-k2+1:m);
%split = sort_ind(~ismember(sort_ind,vert_indices));
split = sort_ind(~ismember(sort_ind,union(vert_indices,singletons)));
split = split(end-k2+1:end);
%split = split(am(split)>0);
split = sort(split,'descend');
%assert(k2 == length(split));
k = max(k2 - num_add,k0); %number of removals 

num_false_add = 0; %number of points which don't get added
y_new = [y_new;zeros(k2,d)];
neighbors_new = [neighbors_new;cell(k2,1)];
% y_new = [y;zeros(k2,d)];
% neighbors_new = [neighbors;cell(k2,1)];

for i=1:k2 
    
    neigh = neighbors_new{split(i)}(1);
  
    t0 = .5;
    new_point = t0*y_new(split(i),:)+(1-t0)*y_new(neigh,:);
    y_new(m+1,:) = new_point;
    if length(neighbors_new{split(i)})>1
        other_neigh = neighbors_new{split(i)}(2);
        neighbors_new{split(i)} = [other_neigh,m+1];
    else
        neighbors_new{split(i)} = [m+1];
    end
    neighbors_new{neigh} = union(setdiff(neighbors_new{neigh},split(i)),m+1);
    neighbors_new{m+1} = [split(i),neigh];

    
%     [net_edges_row_ind, net_edges_ind_col] = ind2sub(size(net_edges),find(net_edges==split(i)));
%     [edges_row_ind, edges_ind_col] = ind2sub(size(net_edges),find(edges==split(i)));
%     assert(length(net_edges_row_ind)==2);
%     for j=1:2
%         if sum(net_edges(net_edges_row_ind(j)) == [split(i),neigh])==2
%             net_edges(net_edges_row_ind(j)) = [split(i),split(i)+1];
%         end
%         if sum(edges(edges_row_ind(j)) == [split(i),neigh])==2
%             edges(edges_row_ind(j)) = [split(i),split(i)+1];
%         end
%     end
     %edges(intersect(edge_inds(split(i)),edge_inds(neigh)),:) = [];
     e_rem_ind = intersect(y_net_edge_inds{split(i)},y_net_edge_inds{neigh});
     net_edges(e_rem_ind,:) = [];  %remove edge between split(i) and neigh
     y_net_edge_inds{split(i)} = setdiff(y_net_edge_inds{split(i)},e_rem_ind);
     y_net_edge_inds{neigh} = setdiff(y_net_edge_inds{neigh},e_rem_ind);
     for j=1:m
        %y_net_edge_inds{j} = y_net_edge_inds{j} - double(y_net_edge_inds{j}>e_rem_ind);
        if ~isempty(y_net_edge_inds{j})
            y_net_edge_inds{j} = y_net_edge_inds{j} - sum(y_net_edge_inds{j}>e_rem_ind',1);  %in degenerate cases may happen that e_rem_ind has length>1
        end
     end
     
     
%     if neigh<split(i)
%         net_edges = [net_edges;[split(i)+1,neigh];[split(i)+1,split(i)]];
%         edges = [edges;[split(i)+1,neigh];[split(i)+1,split(i)]];
%     else
%         net_edges = [net_edges;[neigh,split(i+1)];[split(i)+1,split(i)]];
%         edges = [edges;[neigh,split(i)+1];[split(i)+1,split(i)]];
%     end
    net_edges = [net_edges;[neigh,m+1];[split(i),m+1]]; % add edges from new point to split(i) and neigh
    num_net_e = length(net_edges);
    %edges = [edges;[m+1,neigh];[m+1,split(i)]];
    %num_e = length(edges);
    
    y_net_edge_inds{split(i)} = [y_net_edge_inds{split(i)},num_net_e];
    y_net_edge_inds{neigh} = [y_net_edge_inds{neigh},num_net_e-1];
    y_net_edge_inds = [y_net_edge_inds;{[num_net_e-1,num_net_e]}];
    m = m + 1;
    hold on;
    scatter(new_point(:,1), new_point(:,2),100,'red', 'fill');
end
    

if d==2
    scatter(y_new(removals,1),y_new(removals,2),100,'black','fill');
    % else if d==3
    %         scatter3(y(removals,1),y(removals,2),y(removals,3),'black','fill');
    %     end
end
hold off;

%removals = sort(removals(1:end-num_false_add),'descend');
b_new = []; z_new = [];
%[y,z_new,b_new] = removePoints(y,z_new,b_new,removals);
%k = length(removals);
k = min(k,length(removals));
removals = removals(1:k);
removals = sort(removals,'descend');
for i=1:k
    rem_ind = removals(i);
    y_new(rem_ind,:) = [];
    
    vert_indices = vert_indices - double(vert_indices>=rem_ind);
    
    %edges(edge_inds(rem_ind),:) = [];
    rem_e_inds = y_net_edge_inds{rem_ind};
    if ~isempty(rem_e_inds) %removal has degree 2: account for removing edges
        rem_edges = net_edges(rem_e_inds,:);
        neigh1 = setdiff(rem_edges(1,:),rem_ind);
        neigh2 = setdiff(rem_edges(2,:),rem_ind);
        %neighs = setdiff(union(net_edges(net_edge_inds{rem_ind},:),[]),rem_ind);
        net_edges = [net_edges;[neigh1,neigh2]];
        net_edges(rem_e_inds,:) = [];
        net_edges = net_edges - double(net_edges>=rem_ind);
        
        %y_net_edge_inds{neigh1} = setdiff(y_net_edge_inds{neigh1},rem_e_inds(1));
        %y_net_edge_inds{neigh2} = setdiff(y_net_edge_inds{neigh2},rem_e_inds(2));
        y_net_edge_inds{neigh1} = setdiff(y_net_edge_inds{neigh1},rem_e_inds);
        y_net_edge_inds{neigh2} = setdiff(y_net_edge_inds{neigh2},rem_e_inds);
        for j=1:length(y_net_edge_inds)
            %y_net_edge_inds{j} = y_net_edge_inds{j} - double(y_net_edge_inds{j}>rem_e_inds(1))-double(y_net_edge_inds{j}>rem_e_inds(2));
            if ~isempty(y_net_edge_inds{j})
                y_net_edge_inds{j} = y_net_edge_inds{j} - sum(y_net_edge_inds{j}>rem_e_inds',1);
            end
        end
        y_net_edge_inds{neigh1} = [y_net_edge_inds{neigh1},length(net_edges)];
        y_net_edge_inds{neigh2} = [y_net_edge_inds{neigh2},length(net_edges)];
    else %removal is singleton
        net_edges = net_edges - double(net_edges>=rem_ind);
    end
    y_net_edge_inds = {y_net_edge_inds{1:rem_ind-1},y_net_edge_inds{rem_ind+1:end}};
    
end

% neighbors = cell(m-k,1);
% neighbors(vert_indices) = vert_neighs;
% for i=1:length(arcs)
%     for j=2:length(arcs{i})-1
%         neighbors{arcs{i}(j)} = [arcs{i}(j-1),arcs{i}(j+1)];
%     end
% end
% 
% cut_indices(cut_indices<1)=[];

%new_struct = {cut_indices,vert_indices,neighbors,arcs,net_edges,arc_inds,net_edge_inds};


end

