function [y,net_edges,edges,edge_weights,iter,energy,dists,next] = onut(y0,net_edges0,x,mass,lambda,alpha,alpha_2,theta,...
    tol,rho0,max_m,max_avg_turn,normalize_data,plot_bool,delta,max_leng)
%ONUT computes optimal network for universal transport
%   This function computes approximate local minimizers to the universal
%   transportation functional.
%INPUTS:
% y0 - locations of points on initial network, given as mxd vector
%   (m-number of points, d-dimension)
% net_edges - a list of edges for the initial network. Each edge contains
%   the indices of the points in y0 that form the edge. It is given as a
%   n_e x 2 vector, where n_e is the number of edges
% x - the n data points, given as nxd vector
% mass - the data masses, given as 1xn vector
% lambda - parameter in ONST functional (value > 0)
% alpha - cost of travelling along network, value between 0 and 1.
% tol - the tolerance for convergence of the energy. 
% rho0 - the initial value of rho, a parameter in the ADMM algorithm (it
%   may be adaptively changed based on the residuals)
% max_m - the maximum number of points on the network. This criteria may
%   be dominated by criteria given by max_leng parameter.
% max_avg_turn - not currently implemented to have effect -- but would be the
%   maximum desired average turning angle, in degrees, of the network, and 
%   which serves as a criteria for when the network is well enough resolved 
% normalize_data - value 0 or 1 whether to normalize the data x
% plot_bool - value 0 or 1 whether to plot after every iteration
% delta - longest length of edges that may be added using heuristic, value
%   > 0
% max_leng - maximum length of edges on network, edges longer will be
% further resolved. This criteria dominates the max_m criteria.
%OUTPUTS:
% y - locations of points on final network, given as mxd vector
%   (m-number of points, d-dimension)
% net_edges - a list of edges for the initial network. Each edge contains
%   the indices of the points in y0 that form the edge. It is given as a
%   n_e x 2 vector, where n_e is the number of edges
% edges - a list of edges on the graph (XuY), whose vertices consist of all
%   data points 'x', and all network points 'y'.
% edge_weights - a list of weights corresponding to 'edges' (of the graph
%   XuY). This vector contains the total mass traveling through each edge,
%   'm_{a,b}' as defined in the paper.
% iter - the total number of (outer) iterations of the onut algorithm.
% energy - the discrete energy of the final network (y,net_edges).
% dists - the geodesic (shortest path) distances between all m+n vertices 
%   on the complete graph XuY, with weights equal to length of the edge,
%   times alpha if the edge belongs to the network.
% next - an (m+n)x(m+n) matrix of indices indicating the shortest paths
%   between all pairs of vertices on the graph XuY. the entry at next(i,j)
%   gives the point that should be visited next if going from point i to
%   point j. 

FIX_ENDPTS = 0;
MAX_INNER_ITER = 100;
r = 4; %param for modifySpacing
n = length(x(:,1));

if isempty(max_m)
    max_m = n/2;
end
if normalize_data
    x_mean = mass*x;
    x_var = mass*sum((x-repmat(x_mean,n,1)).^2,2);
    x = (x-repmat(x_mean,n,1))/sqrt(x_var);
    m = length(y0(:,1));
    y0 = (y0 - repmat(x_mean,m,1))/sqrt(x_var);
end

if isempty(tol), tol = 10^-3; end
num_d = ceil(-log(tol)/log(10)); %number of decimals to display of energy
if isempty(rho0)
    rho = lambda;
else
    rho = rho0;
end

y = y0; z = []; b = [];
rgb_c = [0,.8,0]; ls = '-';
energy_prev = inf;
iter = 0; check_top = 0; e_heur = 0; more_prec = 0; n_add = 1;
%fig1 = figure;
%plotNet(x,y,net_edges0,rgb_c,ls);
%pause;
[edges,edge_costs,edge_weights,net_edges,y_net_edge_inds,neighbors,y,z,b,next,top_change,dists] = computeEdgeWeights(y,x,alpha,lambda,net_edges0,mass,check_top,e_heur,delta);
energy = calculateEnergy(y,x,edges,edge_costs);
energy01 = calculateEnergyOld(y,x, net_edges, mass, alpha_2, theta, lambda);
energy_1 = calculateEnergyNew(y,x, net_edges, mass, alpha_2, theta, lambda);
%fig2 = figure;
%plotGraph(x,y,edges,edge_weights,lambda,alpha,20,mass,rgb_c,ls);
%energy0 = energy; %will be energy after ADMM
%fprintf(['\n iter = %d       E = %8.',num2str(num_d),'E'],iter,energy);
energy0 = inf;
while energy_prev-energy>tol*energy || energy0-energy>tol*energy || energy0<energy || energy_prev<energy || ~check_top || top_change || n_add>0
    energy_prev = energy;
    iter = iter + 1;
    inner_iter = 0;
    energy0 = inf; 
    fi = 1;
    y_new = y; z_new = z; b_new = b;
    while (energy0>energy_prev || (more_prec && (energy0-energy_prev>tol*energy0)) || inner_iter<10)...
            && inner_iter<MAX_INNER_ITER 
        z_prev = z_new;
        [y_new,z_new,b_new] = updateyzb(y_new,z_new,b_new,x,rho,edges,edge_costs,fi,FIX_ENDPTS);
        fi = 0;
        energy0 = calculateEnergy(y_new,x,edges,edge_costs);
%         energy01 = calculateEnergyOld(y_new,x, net_edges, mass, alpha_2, theta, lambda);
%         energy_1 = calculateEnergyNew(y_new,x, net_edges, mass, alpha_2, theta, lambda);
%         fprintf('Old energy0: %f, %f, %d\n', energy0, energy01,inner_iter);
        %figure(fig1); plotNet(x,y,net_edges,rgb_c,ls); drawnow;
        %figure(fig2); plotGraph(x,y_new,edges,edge_weights,lambda,alpha,20,mass,rgb_c,ls); drawnow;
        inner_iter = inner_iter + 1;
%         if energy0>energy_prev && inner_iter>10
%             y_new = y; z_new = z; b_new = b;
%             rho = rho*1.3;
%             fprintf('    rho = %3.2E',rho);
%         end
    end
    if energy0>energy_prev %inner_iter==MAX_INNER_ITER &&
        fprintf('\n MAX INNER ITERATIONS (%d) REACHED, ENERGY DID NOT DECREASE',MAX_INNER_ITER);
    else
        y = y_new; z = z_new; b=b_new;
    end
    if inner_iter>1
        %fprintf('\n %d inner ADMM iterations',inner_iter);
        %fprintf('\n %d inner iters',inner_iter);
        if ~isempty(z_new) && length(z_prev(:,1))==length(z_new(:,1))
            [r_pri,r_dual] = computeResiduals(x,y,edges,z_new,z_prev,rho);
            if r_pri>5*r_dual, rho = rho*1.3;
                if ~isempty(rho0), %fprintf('    rho = %3.2E',rho); 
                end
            end
            if 5*r_pri<r_dual, rho = rho/2;
                if ~isempty(rho0), %fprintf('    rho = %3.2E',rho);
                end
            end
        end
    end
    if more_prec == 1
        %fprintf('  (more precision in minimization obtained)');
        more_prec = 0;
    else
        
    end
    
    %fprintf(['\n iter = %d       E = %8.',num2str(num_d),'E'],iter,energy0);
    %if check_top && ~top_change
    [y,net_edges] = resolveLongEdges(y,net_edges,max_leng);
    %end
    %[y,net_edges] = mergeVertices(y,net_edges,delta); %ad-hoc routine: if
    %used change 'delta' parameter name. delta is now used as a parameter for checking edge additions
    
    %if iter>=15 
    if iter>=0 
        check_top = 1;
        if top_change %only test heuristic if no trivial edge changes
            e_heur = 0;
        else
            if iter>1, e_heur = 1; end
        end
    else
        check_top = 0;
    end
    
    [edges,edge_costs,edge_weights,net_edges,y_net_edge_inds,neighbors,y,z,b,next,top_change,dists] = computeEdgeWeights(y,x,alpha,lambda,net_edges,mass,check_top,e_heur,delta);
    energy = calculateEnergy(y,x,edges,edge_costs);
    energy01 = calculateEnergyOld(y,x, net_edges, mass, alpha_2, theta, lambda);
    energy_1 = calculateEnergyNew(y,x, net_edges, mass, alpha_2, theta, lambda);%In each place of calculating E, we add a calculation of E in the new version
    if energy0-energy < 10*tol*energy
        more_prec = 1;
    end

    if mod(iter,5)==0
        if mod(iter,10)==0
            eo = 1;
        else
            eo = 0;
        end
        n_add = numPointsAdd(y,net_edges,y_net_edge_inds,neighbors,max_avg_turn,max_m,max_leng);
        top_struct = {edges,edge_weights,net_edges,y_net_edge_inds,neighbors};
        [y,z,b,net_edges] = modifySpacingAllComps(y,r,top_struct,n_add,max_m,eo,1,FIX_ENDPTS);
        [edges,edge_costs,edge_weights,net_edges,~,neighbors,y,z,b,next,top_change,dists] = computeEdgeWeights(y,x,alpha,lambda,net_edges,mass,check_top,e_heur,delta);
        energy = calculateEnergy(y,x,edges,edge_costs);
        %fprintf(['\n       E = %8.',num2str(num_d),'E      m = %d'],energy,length(y(:,1)));
    end
    if plot_bool
        %figure(fig1); plotNet(x,y,net_edges,rgb_c,ls); drawnow;
        %figure(fig2); plotGraph(x,y,edges,edge_weights,lambda,alpha,20,mass,rgb_c,ls); drawnow;
        %pause;
    end
    %fprintf(['\n iter = %d       E = %8.',num2str(num_d),'E'],iter,energy);
end

if normalize_data %un-normalize data
    x = sqrt(x_var)*x + repmat(x_mean,n,1);
    m = length(y(:,1));
    y = sqrt(x_var)*y + repmat(x_mean,m,1);
end
end



function [r_pri,r_dual] = computeResiduals(x,y_new,edges,z_new,z,rho) %the output of this function is used to modify rho
[m,d] = size(y_new);
num_intern_edges = sum(edges(:,1)<=m);
rel_edge_ind = edges(:,2)<=m;
edges = edges(rel_edge_ind,:);
num_e = length(edges(:,1));
z0_1 = y_new(edges(1:num_intern_edges,2),:) - y_new(edges(1:num_intern_edges,1),:);
z0_2 = y_new(edges(num_intern_edges+1:num_e,2),:); %assumes y indices are second coordinate in edges(i)%
c = [zeros(num_intern_edges,d);-x(edges(num_intern_edges+1:num_e,1)-m,:)];
r_pri = norm([z0_1;z0_2] + c - z_new);
%r_pri = norm(y_new(2:end,:)-y_new(1:end-1,:)-z_new);
z_diff = z_new - z;
%r_dual = rho*norm([z_diff(1,:);z_diff(1:end-1,:)-z_diff(2:end,:);z_diff(end,:)]);
DTzdiff = zeros(m,d);
for i=1:num_e
    if edges(i,1)<= m
        DTzdiff(edges(i,1)) = DTzdiff(edges(i,1)) - z_diff(i);
    end
    DTzdiff(edges(i,2)) = DTzdiff(edges(i,2)) + z_diff(i);
end
r_dual = rho*norm(DTzdiff);
end

