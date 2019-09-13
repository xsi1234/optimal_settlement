function E = calculateEnergyNew(Y,X, net_edges,M, alpha, theta, lambda)
    Adj = zeros(size(Y,1));
    for i = 1:size(net_edges, 1)
        Adj(net_edges(i,1), net_edges(i,2)) = 1;
    end
    Adj = or(Adj, Adj');
    l = size(X,1);
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
                   theta1 = acos(alpha);
                   beta = acos(dot(xi-yk, yj-yk)/(norm(xi-yk)*norm(yj-yk)));
                   gamma = acos(dot(xi-yj, yk-yj)/(norm(xi-yj)*norm(yj-yk)));
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
    [dists, next] = floyd(D_mat);
    dists = dists(ly+1:l+ly, ly+1:ly+l);
    E1 = sum(sum(dists .* repmat(M, l, 1) .* repmat(M', 1, l)));
    E2 = sum(sum(Adj .* squareform(pdist(Y))))/2;
    K_mat = zeros(l);
    for i = 1:l
        for j = 1:l
            if i ~= j
                K_mat(i,j) = computeK(X(i, :) - X(j, :), 0.05);
            end
        end
    end
    E3 = ((M * K_mat) .^ (3-1))*M';
    E = lambda*E2 + E1+ theta*E3;
    fprintf(['E =%.6f   E1 = %.6f   E2 = %.6f   E3 = %.6f\n'],E, lambda*E2, E1, theta*E3);
end    
    
    
function [dists,next] = floyd(weight_mat)
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

function y = computeK(x, sigma)
y = exp(-x*x'/(sigma^2))/sigma;
end