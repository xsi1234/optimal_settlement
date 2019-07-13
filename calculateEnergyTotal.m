function E = calculateEnergyTotal(y, Adj, dist_mat, K_mat, theta, lambda)
    n = size(K_mat,1);
    Adj = triu(Adj); 
    network_cost = Adj .* squareform(pdist(y));
    E1 = sum(network_cost, 'all');
    E2 = sum(dist_mat, 'all') * (1/n)^2;
    E3 = sum(K_mat, 'all') * (1/n)^2;
    E = lambda * E1 + E2 + theta * E3;
end