function E = calculateEnergyTotal(y, Adj, dist_mat,old_dist_mat, K_mat, theta, lambda, p, M)
    total_mass = sum(M);
    n = size(K_mat,1);
    for i = 1:n
        K_mat(i,i) = 0;
    end
    Adj = triu(Adj); 
    network_cost = Adj .* squareform(pdist(y));
    E1 = sum(sum(network_cost));
    E2 = sum(sum(dist_mat .* repmat(M, n, 1) .* repmat(M', 1, n)));
    E2_old = sum(sum(old_dist_mat .* repmat(M, n, 1) .* repmat(M', 1, n)));
    E3 = ((M * K_mat) .^ (p-1))*M';
    E = lambda * E1 + E2 + theta * E3;
    fprintf(['E =%.6f   E1 = %.4f   E2 = %.4f   old_E2 = %.4f   E3 = %.4f \n'],E, lambda*E1 , E2, E2_old, theta*E3);
end