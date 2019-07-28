%Wrapping function for forward Euler Scheme. h = step size, p = kernel
%dimension. Returns the positions of particles after iter_num iterations.
function X_res = euler_wrap(X, Y, Adj, M, iter_num, p, h0, alpha, theta, sigma, lambda, color_mat)
 %   figure();
 %   hold on;
 %   plot_network(Y, Adj);
    tic
    h = h0;
    X_record = [];
    E_last = inf;
    for i = 1:iter_num
        if h < h0*0.001
            break
        end
        [E,grad,X] = euler_iter(X, M, Y, Adj, h, p, theta, alpha,sigma,lambda);
        fprintf(['\n iter = %d   E = %.8f'],i,E);
        if E > E_last
            h = h * 1/2;
        end
        E_last = E;
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

function plot_particles(X, X_record, color_mat)
%     for i = 1:size(X_record, 1)
%         plot(X_record(i,1), X_record(i,2), 'marker','.','color', color_mat(mod(i-1,size(X,1))+1,:));
%     end
    for i = 1:size(X,1)
        plot(X(i,1), X(i,2),'marker','o','markersize',5, 'color',color_mat(i,:));
    end
end

function y = computeK(x, sigma)
y = exp(-x*x'/(sigma^2))/sigma;
end

function y = computeKgrad(x, sigma)
y = -2* x * computeK(x,sigma)/(sigma^2); 
end

% function y = computeK(x, sigma)
% y = exp(-norm(x)/(sigma^2))/sigma;
% if y < 1e-20
%     y=0;
% end
% end
% 
% function y = computeKgrad(x, sigma)
% y = -(x /norm(x))*computeK(x, sigma)/(sigma^2); 
% end

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
        next_i = next(ly+k, ly+i);
        if next_i ~= 0
            flag = 0;
            zi = XY_union(next_i, :);
            next_2 = next(next_i, ly+i);
            if next_2 > 0 && next_2 <= ly && next_i <=ly
                z2 = XY_union(next_2, :);
                v1 = xk - zi;
                v2 = zi - z2;
                beta = acosd(dot(v1,v2)/(norm(v1)*norm(v2)));
                theta = acosd(alpha);
                if beta > theta
                    v = -(z2-zi)/norm(z2-zi);
                    xp = dot(xk-zi,v)*v+zi; 
                    wi = xp + v * min(norm(xk - xp) * cot(theta), norm(z2 - xp));
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
function [E,grad,X_new] = euler_iter(X, M, Y, Adj, h, p, theta, alpha,sigma, lambda)
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
    [dists, next] = floyd1(D_mat);
    for i = 1:l
        %theta*nextstep_E2(X, l, i, M, p, K_mat_1,sigma)
        X_new(i,:) = X(i,:)+h*(nextstep_E1(l, ly, XY_union, next, i, M, alpha)+theta*nextstep_E2(X, l, i, M, p, K_mat_1,sigma));
    end
    grad = X_new - X;
    E = calculateEnergyTotal(Y, Adj, dists(ly+1:l+ly, ly+1:ly+l), K_mat, theta, lambda, p, M);
    %grad = grad / norm(grad);
    %step_len = sum(vecnorm(X_new-X))/h;
    %step_len
end