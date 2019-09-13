function E = calculateEnergyTotal(Y, Adj, X, K_mat, theta, lambda, p, M,alpha)
    total_mass = sum(M);
    n = size(K_mat,1);
    for i = 1:n
        K_mat(i,i) = 0;
    end
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
    [dists, ~] = floyd(D_mat);
    [old_dists, ~] = floyd(old_D_mat);
    dist_mat = dists(ly+1:l+ly, ly+1:ly+l);
    old_dist_mat = old_dists(ly+1:l+ly, ly+1:ly+l);
    Adj = triu(Adj); 
    network_cost = Adj .* squareform(pdist(Y));
    E1 = sum(sum(network_cost));
    E2 = sum(sum(dist_mat .* repmat(M, n, 1) .* repmat(M', 1, n)));
    E2_old = sum(sum(old_dist_mat .* repmat(M, n, 1) .* repmat(M', 1, n)));
    E3 = ((M * K_mat) .^ (p-1))*M';
    E = lambda * E1 + E2 + theta * E3;
    fprintf(['E =%.6f   E1 = %.6f   E2 = %.6f   E3 = %.6f   old_E2 = %.6f   old_E = %.6f\n'],E, lambda*E1 , E2, theta*E3, E2_old, lambda * E1 + E2_old + theta * E3);
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