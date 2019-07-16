function E = calculateEnergyTotal(y, Adj, dist_mat, K_mat, theta, lambda, p, M)
    total_mass = sum(M);
    n = size(K_mat,1);
    Adj = triu(Adj); 
    network_cost = Adj .* squareform(pdist(y));
    E1 = sum(sum(network_cost));
    E2 = sum(sum(dist_mat)) * (total_mass/n)^2;
    E3 = sum(sum(K_mat).^2) * (total_mass/n)^p;
    fprintf(['  E1 = %.4f   E2 = %.4f    E3 = %.4f'],lambda*E1,E2, theta*E3);
    E = lambda * E1 + E2 + theta * E3;
end