function p = floyd(D)
N = length(D(1,:));
dists = D;
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
next1 = repmat(1:N,N,1) - diag(1:N);
for k=1:N
	i2k = repmat(D(:,k), 1, N);
	k2j = repmat(D(k,:), N, 1);
	D1 = min(D, i2k+k2j);
    diff = double(D1-D~=0);
    same = ones(N) - diff;
    next_ik = repmat(next(:,k), 1, N);
    next1 = next1 .* same + diff .* next_ik; 
    D = D1;
end
k = next - next1
p=dists-D;
end