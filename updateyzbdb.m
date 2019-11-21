function [ y_new,z_new,b_new] = updateyzb(y,z,b,x,mass,lambda,rho,neighbors,I,edges,edge_weights,first_iter,FIX_ENDPTS)
%ADMMSTEP
% This function performs one step of the ADMM/split bregman algorithm, updating
% each of y,z,b once. It also takes into account seperate components of y
% given by cut_indices. If multiple components exist, the ADMM step is run
% once on each component with its corresonding data and mass.
% If a particular component has no mass (no corresponding data), then it is
% deleted entirely from y, and z,b are adjusted accordingly.

% INPUTS:
% y-current points for curve
% z-variable for split bregman
% b-variable for split bregman
% x-the data points
% mass-the mass of the data points
% lambda-the coefficient of the length penalty
% rho-coefficient of the penalty for the discrepancy between phi(y) and z
% cut_indices-the last y indices of each component
% neighbors - cell containing the indices of the neighbors of each point in
% y. neighbors{i} is the vector of indices of neighbors of y(i).
% I - if it is nonempty, this one will be used (new projections won't be
% computed) (not used)
% edges - matrix with two columns, each representing indices of oriented
% edges. The first entry of row i represents the starting vertex of edge i,
% and we assume edges(i,1)>=edges(i,2) for all i.
% num_intern_edges - the number of internal edges (edges between points on the
% network). Edge i for i> intern_edges corresponds to an edge between the
% data and the network.
% edge weights - vector with entry i corresponding to weight of edge i.
% first_iter - a boolean indicating whether this is the first iteration of
% ADMM after computing edge weights.
% FIX_ENDPTS - if 1 then updates are done such that endpoints remain fixed.
% Otherwise usual updates are performed.
%
% OUTPUTS:
% I-the assignments of data to closest points in y
% energy-the energy functional evaluated at y_new, not including the
% lambda2 term
% the rest are updated versions of the corresponding inputs
%
% Differs from public version in that remZeroMassPoints takes in FIX_ENDPTS
% parameter. Also has a couple of lines asserting that.

[m,d] = size(y);
num_intern_edges = sum(edges(:,1)<=m); % edges(:,1) <= edges(:,2) elementwise ALSO ASSUMES edges(:,1) is sorted
%z0 = y(2:end,:)-y(1:end-1,:);
% if isempty(z)
%     num_e = 0;
%     for i=1:m
%         num_e = num_e + length(neighbors{i});
%     end
%     num_e = num_e/2;
% else
%     num_e = length(z(:,1));
% end
% z0 = zeros(num_e,d);
% last_ind = 0;
% e_ind = zeros(m,m);
% for i=1:m
%     for j=1:length(neighbors{i})
%         if neighbors{i}(j)>i
%             last_ind = last_ind + 1;
%             z0(last_ind,:) = y(neighbors{i}(j),:) - y(i,:);
%             e_ind(i,neighbors{i}(j)) = last_ind;
%         end
%     end
% end
assert(sum(edges(:,2)>=edges(:,1))==0);
rel_edge_ind = edges(:,2)<=m;
edges = edges(rel_edge_ind,:); %disregard edges between data points
edge_weights = edge_weights(rel_edge_ind,:);
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
if isempty(z)
    prez = z0+b;
    normsprez = sqrt(sum(prez.^2,2)); % computes the euclidean norm of each row vector
    normsprez(normsprez==0)=1;
    normalprez = bsxfun(@rdivide,prez, normsprez);
    z = bsxfun(@times,normalprez, max(normsprez-edge_weights/rho,0));
end
b = b + z0 - z;



%Step 3: update y's
% if isempty(I)
%     [I,~] = dsearchn(y,x);
%     I = I';
% end

%y_ends = [y(1,:);y(end,:)];
%[y_new,z_new,b_new,I,y_mass,cut_indices] = remZeroMassPoints(y,z,b,mass,I,cut_indices,FIX_ENDPTS);
z_new = z; b_new = b;
y_mass = zeros(m,1);
% for i=1:m
%     y_mass(i) = sum(mass(I==i));
% end
% if FIX_ENDPTS
%     assert(isequal(y(1,:),y_ends(1,:)) && isequal(y(end,:),y_ends(2,:)))
%     %if ~isequal(y(1,:),y_ends(1,:))
%     %    y = [y_ends(1,:);y]; end
%     %if ~isequal(y(end,:),y_ends(2,:))
%     %    y = [y;y_ends(2,:)]; end
% end
% [m,~] = size(y);

y_new = updateY(y,z_new,b_new,rho,c,edges,num_intern_edges,FIX_ENDPTS);

% if isempty(cut_indices)
%     y_new = updateY(y,z_new,b_new,rho,c,edges,num_intern_edges,FIX_ENDPTS);
%     
% else
%     y_new = zeros(size(y));
%     cut_indices = union(cut_indices(1),cut_indices); %to avoid repeated values: (why) is it needed?
%     cut_partition = [0;cut_indices;m];
%     
%     %energy1 = 0;
%     
%     for j=2:length(cut_indices)+2   %perform updates on each component. This is parallelizable
%         
%         comp_y_indices = cut_partition(j-1)+1:cut_partition(j);
%         comp_y = y(comp_y_indices,:);
%         %comp_ymass = y_mass(comp_y_indices);
%         %comp_neigh = neighbors(comp_y_indices);
%         %comp_e_ind = e_ind(comp_y_indices,comp_y_indices) - e_ind(comp_y_indices(1),comp_y_indices(2)) + 1;
%         
%         comp_edge_indices = edges(:,2)>=cut_partition(j-1)+1 & edges(:,2)<=cut_partition(j);
%         
%         comp_edges = edges(comp_edge_indices,:) - cut_partition(j-1);
%         comp_num_intern_edges = sum(comp_edges(:,1)<=length(comp_y_indices));
%         
%         comp_c = [zeros(comp_num_intern_edges,d);-x(comp_edges(comp_num_intern_edges+1:end,1)-m+cut_partition(j-1),:)];
%         
%         %comp_xindices = find((I>cut_partition(j-1) & I<=cut_partition(j)));
%         %comp_x = x(comp_xindices,:);
%         %comp_mass = mass(comp_xindices);
%         
%         %comp_I = I(comp_xindices)-cut_partition(j-1)*ones(1,length(comp_xindices));
%         %comp_z_indices = cut_partition(j-1)+1:cut_partition(j)-1;
%         comp_z = z_new(comp_edge_indices,:);
%         comp_b = b_new(comp_edge_indices,:);
%         
%         %comp_y = updateY(comp_y,comp_z,comp_b,comp_x,...
%         %    comp_I,comp_mass,comp_ymass,rho,comp_neigh,comp_e_ind,cut_partition(j-1),FIX_ENDPTS);
%         comp_y = updateY(comp_y,comp_z,comp_b,rho,comp_c,comp_edges,comp_num_intern_edges,FIX_ENDPTS);
%         
%         y_new(comp_y_indices,:) = comp_y;
%         
%         %if length(comp_y(:,1)) > 1
%         %    energy1 = energy1 + lambda*sum(sqrt(sum((comp_y(2:end,:)-comp_y(1:end-1,:)).^2,2)));
%         %end
%         
%     end
%     %energy = energy1 + mass*sum((x-y_new(I,:)).^2,2);
% end
end




function y_new = updateY(y,z,b,rho,c,edges,intern_edges,FIX_ENDPTS)
%Helper function used for updating y on a single component.

[m,d] = size(y);

% x_bar = zeros(m,d);
% for i=1:m
%     %assert(y_mass(i)~=0);
%     x_ind = I==i;
%     %x_bar(i,:) = sum(bsxfun(@times,mass(x_ind),x(x_ind,:)'),2);
%     x_bar(i,:) = mass(x_ind)*x(x_ind,:);
% end
%
% if m==1
%     y_new = bsxfun(@rdivide,x_bar,y_mass);
% else
%     D = 2*y_mass+2*rho*ones(m,1);
%     D(1) = 2*y_mass(1) + rho;
%     D(m) = 2*y_mass(m) + rho;
%     zb1 = [zeros(1,d);z-b];
%     zb2 = [z-b;zeros(1,d)];
%     zb = zb1 - zb2;
%     RHS = 2*x_bar + rho*zb;

%     A = zeros(m,m);
%     B = 2*x_bar;
%     for i=1:m
%         for j=1:length(neighbors{i})
%             neigh = neighbors{i}(j) - offset;
%             A(i,i) = A(i,i) + 1;
%             A(i,neigh) = - 1;
%             if i<neigh
%                 B(i,:) = B(i,:) - rho*(z(e_ind(i,neigh),:) - b(e_ind(i,neigh),:));
%             else
%                 B(i,:) = B(i,:) + rho*(z(e_ind(neigh,i),:) - b(e_ind(neigh,i),:));
%             end
%         end
%     end
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
    %y_new = tridiag((-rho)*ones(m-1,1),D,(-rho)*ones(m-1,1),RHS);
    %         U = zeros(m-1,m);
    %         U(1:m-1,2:m) = eye(m-1);
    %         phi = -eye(m-1,m)+U;
    %         L = zeros(m);
    %         L(2:m,1:m-1) = -rho*eye(m-1);
    %         U = zeros(m);
    %         U(1:m-1,2:m) = -rho*eye(m-1);
    %         A = U + L + diag(2*y_mass+2*rho*ones(m,1));
    %         A(1,1) = 2*y_mass(1) + rho;
    %         A(m,m) = 2*y_mass(m) + rho;
    %         B = 2*x_bar + rho*phi'*(z-b);
    y_new = linsolve(A,B);
else % If endpoints are fixed:
    if m == 2
        y_new = y;
    else
        y_new(2:m-1,:) = tridiag((-rho)*ones(m-3,1),D(2:m-1),(-rho)*ones(m-3,1),...
            RHS(2:m-1,:) + rho*[y(1,:);zeros(m-4,d);y(m,:)]);
        y_new([1,m],:) = y([1,m],:);
    end
end

end



