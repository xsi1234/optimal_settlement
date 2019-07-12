function [y,x,net_edges,edges,edge_weights,iter,energy,next,lambda1,alpha] = mppntest1( y0,cut_indices,net_edges,vert_indices,vert_neighs,arcs,...
    x,mass,lambda1,alpha,tol,rho0,max_m,max_avg_turn,normalize_data,pause_bool,delta,max_leng)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% vert_indices - a vector containing the indices of all points in y with
% order not equal to two (that is, all endpoints and 'intersections').
% These are vertices of the topolgical graph representation of Sigma.
% vert_neighs - a cell containing the indices of all the neighbors of each
% of vertex (in vert_indices).
FIX_ENDPTS = 0;
r = 3; %param for modifySpacing
if isempty(max_m)
    max_m = length(x(:,1))/2;
end
if normalize_data 
    n = length(x(:,1));
    x_mean = mass*x;
    x_var = mass*sum((x-repmat(x_mean,n,1)).^2,2);
    x = (x-repmat(x_mean,n,1))/sqrt(x_var);
    m = length(y0(:,1));
    y0 = (y0 - repmat(x_mean,m,1))/sqrt(x_var);
end

[m,d] = size(y0);
% num_neighs = zeros(m);
% for i=1:m
%     num_neighs = length(neighbors{i});
% end
% vert_indices = find(num_neighs ~= 2);
% all_neighs_coor = [];
% for i=1:length(vert_neighs)
%     all_neighs_coor = [all_neighs_coor,vert_indices(i)*ones(1,length(vert_neighs{i}))];
% end % what if there are singletons?
% [all_vert_neighs,sort_ind] = sort(cell2mat(vert_neighs)); % sorted list of
% % neighbors of vertices. These are the (second to last) endpoints of each
% % arc (edge in topological graph)
% all_neighs_coor = all_neighs_coor(sort_ind); % list with all vertices with
% % multiplicity according to number of neighbors, and order corresponding to
% % all_vert_neighs
%
% assert(mod(length(all_vert_neighs),2) == 0);
% num_arcs = length(all_vert_neighs)/2;
% arcs = cell(num_arcs,1); % cell containing the ordered indices corresponding to each arc (edge in the top graph)
% neighbors = cell(m,1); % cell containing the indices corresponding to the neighbors of each point
% neighbors(vert_indices) = vert_neighs;
%
% for i=1:num_arcs
%     j = all_vert_neighs(2*i-1);
%     prev_neigh = all_neighs_coor(2*i-1);
%     end_ind = all_neighs_coor(2*i);
%     %assert(end_ind == all_vert_neighs(2*i)+1); %check
%     neighbors{j} = [prev_neigh,j+1];
%     arcs{i} = [prev_neigh,j:end_ind];
%     for j=all_vert_neighs(2*i-1)+1:all_vert_neighs(2*i)
%         neighbors{j} = [j-1,j+1];
%     end
% end

if isempty(net_edges)
    
    neighbors = cell(m,1); % cell containing the indices corresponding to the neighbors of each point
    neighbors(vert_indices) = vert_neighs; %elements of neighbors need not be
    % ordered, since net_edges need not be oriented in any particular direction.
    arc_inds = cell(m,1); %entry i gives the arc indices that y_i is contained in
    y_net_edge_inds = cell(m,1); %entry i gives the net_edge indices that y_i neighbors
    y_verts = cell(m,1); %entry i gives the topological vert indices that y_i 'neighbors'
    
    for i=1:length(arcs)
        arc_inds{arcs{i}(1)} = [arc_inds{arcs{i}(1)},i];
        y_verts{arcs{i}(1)} = [y_verts{arcs{i}(1)},arcs{i}(end)];
        for j=2:length(arcs{i})-1
            neighbors{arcs{i}(j)} = [arcs{i}(j-1),arcs{i}(j+1)];
            arc_inds{arcs{i}(j)} = [arc_inds{arcs{i}(j)},i];
            y_verts{arcs{i}(j)} = [y_verts{arcs{i}(j)},setdiff([arcs{i}(1),arcs{i}(end)],arcs{i}(j))];
        end
        arc_inds{arcs{i}(end)} = [arc_inds{arcs{i}(end)},i];
        y_verts{arcs{i}(end)} = [y_verts{arcs{i}(end)},arcs{i}(1)];
    end
    
    num_net_edges = m-length(vert_indices) + length(cell2mat(vert_neighs))/2.0;
    net_edges = zeros(num_net_edges,2);
    k = 1;
    for i=1:m
        for j=1:length(neighbors{i})
            if neighbors{i}(j)>i
                net_edges(k,:) = [i,neighbors{i}(j)];
                y_net_edge_inds{i} = [y_net_edge_inds{i},k];
                y_net_edge_inds{neighbors{i}(j)} = [y_net_edge_inds{neighbors{i}(j)},k];
                k = k + 1;
            end
        end
    end
end

if isempty(tol), tol = 10^-3; end
num_d = ceil(-log(tol)/log(10)); %number of decimals to display of energy
if isempty(rho0)
    rho = lambda1;
else
    rho = rho0;
end

y = y0; z = []; b = [];
rgb_c = [0,.8,0]; ls = '-';
%fig1 = figure;
%plotADP(x,y0,arcs,mass,[],0,rgb_c,ls);
%plotNet(x,y,net_edges,rgb_c,ls);
energy_prev = inf;
iter = 0; check_top = 0;
[edges,edge_costs,edge_weights,net_edges,y_net_edge_inds,neighbors,y,z,b,next] = computeEdgeWeights1(y,x,alpha,lambda1,net_edges,mass,check_top);
%fig2 = figure;
%plotGraph(x,y,arcs,edges,edge_weights,lambda1,alpha,20,mass,rgb_c,ls);
energy = calculateEnergy1(y,x,edges,edge_costs);
%fprintf(['\n iter = %d       E = %8.',num2str(num_d),'E      cont_E = %8.',num2str(num_d),'E'],iter,energy,cont_energy);
fprintf(['\n iter = %d       E = %8.',num2str(num_d),'E'],iter,energy);

while energy_prev-energy>tol*energy || energy_prev<energy || ~check_top || n_add>0
    energy_prev = energy;
    iter = iter + 1;
    if iter==24
    end
    inner_iter = 0;
    energy0 = inf;
    fi = 1;
    while (energy0>energy_prev || energy0 == inf) || inner_iter<10
        z_old = z;
        [y,z,b,cut_indices] = updateyzb1(y,z,b,x,mass,lambda1,rho,cut_indices,neighbors,[],edges,edge_costs,fi,FIX_ENDPTS);
        fi = 0;
        energy0 = calculateEnergy1(y,x,edges,edge_costs);
        %plotADP(x,y,arcs,mass,[],0,rgb_c,ls); drawnow;
        %plotGraph(x,y,arcs,edges,edge_weights,lambda1,alpha,20,mass,rgb_c,ls); drawnow;
        inner_iter = inner_iter + 1;
    end
    if inner_iter>1
%         if ~isempty(z) && length(z_old(:,1))==length(z(:,1))
%             [r_pri,r_dual] = computeResiduals(y,z,z_old,rho);
%             if r_pri>5*r_dual, rho = rho*1.3;
%                 if ~isempty(rho0), fprintf('    rho = %3.2E',rho); end
%             end
%             if 5*r_pri<r_dual, rho = rho/2;
%                 if ~isempty(rho0), fprintf('    rho = %3.2E',rho); end
%             end
%         end
        %fprintf('\n %d inner ADMM iterations',inner_iter);
        %fprintf('\n %d inner iters',inner_iter);
    end
    fprintf(['\n iter = %d       E = %8.',num2str(num_d),'E'],iter,energy0);
    [y,net_edges] = resolveLongEdges(y,net_edges,max_leng);
    [y,net_edges] = mergeVertices(y,net_edges,delta);

    [edges,edge_costs,edge_weights,net_edges,y_net_edge_inds,neighbors,y,z,b,next] = computeEdgeWeights1(y,x,alpha,lambda1,net_edges,mass,check_top);
    energy = calculateEnergy1(y,x,edges,edge_costs);
    %if energy0-energy<tol*10
    if iter>=15 %&& mod(iter,3)==0
        check_top = 1;
    else
        check_top = 0;
    end
    if mod(iter,5)==0
        if mod(iter,10)==0
            eo = 1;
        else
            eo = 0;
        end
        %plotADP(x,y,arcs,mass,[],0,rgb_c,ls); drawnow;
        n_add = numPointsAdd(y,max_avg_turn,max_m,cut_indices);
        %top_struc = {cut_indices,vert_indices,vert_neighs,neighbors,arcs,net_edges,edges,edge_costs,arc_inds,y_net_edge_inds};
        %[y,z,b,top_struc] = modifySpacingAllComps(y,z,b,r,lambda1,alpha,...
        %    top_struc,n_add,eo,1,FIX_ENDPTS);
        %[cut_indices,vert_indices,neighbors,arcs,net_edges,arc_inds,y_net_edge_inds] = top_struc{:};
        top_struct = {edges,edge_weights,net_edges,y_net_edge_inds,neighbors};
        [y,z,b,net_edges] = modifySpacingAllComps1(y,r,top_struct,n_add,eo,1,FIX_ENDPTS);
        %plotADP(x,y,arcs,mass,[],0,rgb_c,ls); drawnow;
        [edges,edge_costs,edge_weights,net_edges,~,neighbors,y,z,b,next] = computeEdgeWeights1(y,x,alpha,lambda1,net_edges,mass,check_top);
        energy = calculateEnergy1(y,x,edges,edge_costs);
        %fprintf(['\n       E = %8.',num2str(num_d),'E      m = %d'],energy,length(y(:,1)));
    end
    %figure(fig1); plotADP(x,y,arcs,mass,[],0,rgb_c,ls); drawnow;
    %figure(fig1); plotNet(x,y,net_edges,rgb_c,ls); drawnow;
    %figure(fig2); plotGraph(x,y,arcs,edges,edge_weights,lambda1,alpha,20,mass,rgb_c,ls); drawnow;
    %pause;
    %fprintf(['\n iter = %d       E = %8.',num2str(num_d),'E'],iter,energy);
end

if normalize_data %un-normalize data
    x = sqrt(x_var)*x + repmat(x_mean,n,1);
    m = length(y(:,1));
    y = sqrt(x_var)*y + repmat(x_mean,m,1);
end

function [r_pri,r_dual] = computeResiduals(y_new,z_new,z,rho)
r_pri = norm(y_new(2:end,:)-y_new(1:end-1,:)-z_new);
z_diff = z_new - z;
r_dual = rho*norm([z_diff(1,:);z_diff(1:end-1,:)-z_diff(2:end,:);z_diff(end,:)]);
end

end

