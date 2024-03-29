function E = calculateEnergyOld(Y,X, net_edges,M, alpha, theta, lambda)
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
    [dists, next] = floyd(D_mat);
    dists = dists(ly+1:l+ly, ly+1:ly+l);
    network_cost = Adj .* squareform(pdist(Y));
    E1 = sum(sum(network_cost))/2;
    E2 = sum(sum(dists .* repmat(M, l, 1) .* repmat(M', 1, l)));
    E = lambda * E1 + E2;
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