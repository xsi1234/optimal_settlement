%Wrapping function for forward Euler Scheme. h = step size, p = kernel
%dimension. Returns the positions of particles after iter_num iterations.
function X_res = euler_wrap(X, Y, Adj, M, iter_num, p, h0, alpha, theta, sigma, lambda, color_mat, net_edges)
 %   figure();
 %   hold on;
 %   plot_network(Y, Adj);
    tic
    h = h0;
    X_record = [];
    l = size(X,1)
    K_mat = zeros();
    for i = 1:l
        for j = 1:l
            K_mat(i,j) = computeK(X(i, :) - X(j, :), sigma);
        end
    end
    E_last = calculateEnergyTotal(Y, Adj, X, K_mat, theta, lambda, p, M,alpha);
    for i = 1:iter_num
        if h < 0.00001
            plotting1 = X+grad1/30;
            plotting2 = X+grad2/30;
            plotting_sum = X+(grad1+grad2)/30;
            plot_particles(X, color_mat);
            for i = 1:l
                plot([X(i,1) plotting1(i,1)],[X(i,2) plotting1(i,2)], 'LineWidth', 1, 'color', 'b');
                plot([X(i,1) plotting2(i,1)],[X(i,2) plotting2(i,2)], 'LineWidth', 1, 'color', 'g');
                plot([X(i,1) plotting_sum(i,1)],[X(i,2) plotting_sum(i,2)], 'LineWidth', 1, 'color', 'k');
            end
            drawnow;
            break;
        end
        [E,grad1, grad2,X_new] = euler_iter(X, M, Y, Adj, h, p, theta, alpha,sigma,lambda);
        if E > E_last
            h = h * 1/2;
            continue;
        end
        %fprintf(['E =%.6f   E_last = %.6f\n'],E, E_last);
        fprintf("Accepted new configuration of X \n");
        X = X_new;
        E_last = E;
        h = h0;
        if mod(i, ceil(iter_num/20)) == 0%We allow at most 20 traces per point to appear on the plot, making it less messy
            X_record = [X_record', X']';
        end
    end
    toc
    X_res = X;
%    plot_particles(X, X_record, color_mat);
%    plot_grad(X, grad, color_mat);
%    drawnow;
end
%Plotting functions
function plot_grad(X, grad, color_mat)
    %grad = normr(grad);
    K = X+grad;
    for i = 1:size(X,1)
        plot([X(i,1) K(i,1)], [X(i,2) K(i,2)], 'color', color_mat(i,:));
    end
end
function plot_network (Y, Adj)
%     for i = 1:size(X,1)
%         plot(X(i,1), X(i,2),'marker','o','markersize', 8,'color',color_mat(i,:));
%     end
    plot(Y(:,1), Y(:,2), '+');
    for i = 1:size(Y,1)
        for j = 1:size(Y,1)
             if Adj(i,j) ~= 0
                 plot([Y(i,1) Y(j,1)],[Y(i,2) Y(j,2)], 'LineWidth', 3, 'color', 'k');
             end
         end
    end
end

function plot_particles(X, color_mat)
%     for i = 1:size(X_record, 1)
%         plot(X_record(i,1), X_record(i,2), 'marker','.','color', color_mat(mod(i-1,size(X,1))+1,:));
%     end
    for i = 1:size(X,1)
        plot(X(i,1), X(i,2),'marker','o','markersize',5, 'color',color_mat(i,:));
    end
end

% function y = computeK(x, sigma)
% y = exp(-x*x'/(sigma^2))/sigma^2;
% end
% 
% function y = computeKgrad(x, sigma)
% y = -2* x * computeK(x,sigma)/(sigma^2); 
% end

%  function y = computeK(x, sigma)
%  y = exp(-norm(x)/(sigma))/sigma;
%  end
% 
%  function y = computeKgrad(x, sigma)
%  if  norm(x) < 1e-10
%      y=[0,0];
%  else
%    y = -(x /norm(x))*computeK(x, sigma)/(sigma);
%  end;
%  end
% 
  function y = computeK(x, sigma)
  if  (norm(x) > sigma) | (norm(x) < 1e-8)
    y = 0;
  else
    y=-1/(sigma^2)*log(norm(x)/sigma);
  end
  end

  function y = computeKgrad(x, sigma)
  if  (norm(x) > sigma) | (norm(x) < 1e-8)
      y=[0,0];
  else
    y = -(x /(x*x'))/(sigma^2);
  end
  end

%    function y = computeK(x, sigma)
%  if  norm(x) < 1e-10
%  y = 0;
%  else
%  y=1/(x*x');
%  end;
%  end
%
%  function y = computeKgrad(x, sigma)
%  if  norm(x) < 1e-7
%      y=[0,0];
%  else
%    y = -x /((x*x')^2);
%  end;
%  end

%Computing the convolving term for a given index k
function U=nextstep_E2(X, l, k, M, p, K_mat,sigma)
    K_grad_kth = zeros(l, 2);
    for i = 1:l
        if i~=k
            K_grad_kth(i, :) = computeKgrad(X(k, :)-X(i, :),sigma);
        end
    end
    u_component_1 = M * K_grad_kth .* (K_mat(k));
    M_1 = M .* K_mat;
    u_component_2 = M_1 * K_grad_kth;
    U = -(p-1)*(u_component_1+u_component_2);
end

%Floyd-Warshall(Not vectorized version)
function [dists,next] = floyd1(weight_mat)
N = length(weight_mat(1,:));
dists = weight_mat;
next = repmat(1:N,N,1) - diag(1:N);

for k=1:N
    for i=1:N
        for j=1:N
            if dists(i,k) + dists(k,j) < dists(i,j)
                dists(i,j) = dists(i,k) + dists(k,j);
                next(i,j) = next(i,k);
            end
        end
    end
end
end
%Floyd-Warshall (Vectorized version)
function [D,next] = floyd(weight_mat)
N = length(weight_mat(1,:));
D = weight_mat;
next = repmat(1:N,N,1) - diag(1:N);
for k=1:N
	i2k = repmat(D(:,k), 1, N);
 	k2j = repmat(D(k,:), N, 1);
 	D1 = min(D, i2k+k2j);
    diff = double(D1-D~=0);
    same = ones(N) - diff;
    next_ik = repmat(next(:,k), 1, N);
    next = next .* same + diff .* next_ik; 
    D = D1;
end
end

%Computing the gradient of D term for a given index k
function V = nextstep_E1(lx, ly, XY_union, next, k, M, alpha)
    D_grad_kth = zeros(lx, 2);
    xk = XY_union(ly+k,:);
    for i = 1:lx
        next_i = next(ly+k, ly+i);%index of z1
        if next_i ~= 0
            flag = 0;
            zi = XY_union(next_i, :);%access z1
            next_2 = next(next_i, ly+i);%index of z2
            if next_2 > 0 && next_2 <= ly && next_i <=ly
                z2 = XY_union(next_2, :);%access z2
                v1 = xk - zi;
                v2 = zi - z2;
                beta = acos(dot(-v1,v2)/(norm(v1)*norm(v2)));%compute angle between v1 and v2¡¢
                gamma = acos(dot(xk-z2, v2)/norm(xk-z2)*norm(v2));
                theta = acos(alpha);
                if beta < theta && gamma < pi-theta
                    v = -(z2-zi)/norm(z2-zi);%normalized vector from z2 to zi
                    xp = dot(xk-zi,v)*v+zi; %projection of xk onto the vector;
                    wi = xp + v * norm(xk - xp) * cot(theta);%compute the "should be there? station
                    D_grad_kth(i,:) = (xk - wi) / norm(wi - xk);
                    flag = 1;
                end
            end
            if flag == 0
                D_grad_kth(i, :) = (xk-zi)/norm(xk-zi);
            end
        end
    end
    V = -2 * M * D_grad_kth;
end
    
% function V = nextstep_E1(lx, ly, XY_union, next, k, M, alpha)
%     D_grad_kth = zeros(lx, 2);
%     xk = XY_union(ly+k,:);
%     for i = 1:lx
%         next_i = next(ly+k, ly+i);
%         if next_i ~= 0
%             zi = XY_union(next_i, :);
%             D_grad_kth(i, :) = (xk-zi)/norm(xk-zi);
%         end
%     end
%     V = -2 * M * D_grad_kth;
% end


%One iteration of forward Euler Scheme
function [E,grad_E1,grad_E2,X_new] = euler_iter(X, M, Y, Adj, h, p, theta, alpha,sigma, lambda)
    l = size(X,1);
    X_new = zeros(size(X));
%     K_mat = arrayfun(@(i) computeK(X(mod(i,l)+1,:)-X(ceil(i/l),:), sigma), 1:l*l);
%     K_mat = reshape(K_mat, [l,l]);
%     K_mat
    K_mat = zeros(l);
    for i = 1:l
        for j = 1:l
            K_mat(i,j) = computeK(X(i, :) - X(j, :), sigma);
        end
    end
    K_mat_1 = (M * K_mat).^(p-2);
    ly = size(Y, 1);
    XY_union = cat(1, Y, X);
    D_mat = squareform(pdist(XY_union));
    mask = Adj*(1-alpha);
    mask(l+ly, l+ly) = 0;
    D_mat = D_mat - mask .* D_mat;
    old_D_mat = D_mat;
    for k = 1:ly
        for j = 1:ly
            if Adj(k,j) == 1
                yk = Y(k,:);
                yj = Y(j,:);
                for i= 1:l
                   xi = X(i, :);
                   beta = acos(dot(xi-yk, yj-yk)/(norm(xi-yk)*norm(yj-yk)));
                   gamma = acos(dot(xi-yj, yk-yj)/(norm(xi-yj)*norm(yj-yk)));
                   theta1 = acos(alpha);
                   if beta < theta1 && gamma < pi - theta1
                        v1 = xi - yk;
                        v2 = yk - yj;
                        v = -(yj-yk)/norm(yj-yk);%normalized vector from z2 to zi
                        xp = dot(xi-yk,v)*v+yk; %projection of xk onto the vector;
                        px = xp + v * norm(xi - xp) * cot(theta);
                        d = norm(xi-px) + alpha * norm(px-yk);
                        D_mat(ly+i, k) = d;
                   end
                end    
            end
        end
    end
    [dists, next] = floyd1(D_mat);
    [old_dists, old_next] = floyd1(old_D_mat);
    grad_E1 = zeros(l, 2);
    grad_E2 = zeros(l, 2);
    for i = 1:l
        E1_force = nextstep_E1(l, ly, XY_union, next, i, M, alpha);
        E2_force = theta*nextstep_E2(X, l, i, M, p, K_mat_1,sigma);
        grad_E1(i,:) = E1_force;
        grad_E2(i,:) = E2_force;
        X_new(i,:) = X(i,:)+h*E1_force+h*E2_force;
    end
    new_K_mat = zeros(l);
    for i = 1:l
        for j = 1:l
            new_K_mat(i,j) = computeK(X_new(i, :) - X_new(j, :), sigma);
        end
    end
    E = calculateEnergyTotal(Y, Adj, X_new, new_K_mat, theta, lambda, p, M,alpha);
    %grad = grad / norm(grad);
    %step_len = sum(vecnorm(X_new-X))/h;
    %step_len
end
