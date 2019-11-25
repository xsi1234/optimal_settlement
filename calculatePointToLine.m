function D = calculatePointToLine(X, Y, Adj)
    edge = [];
    for i = 1:size(Adj,1)
        for j = 1:i
            if Adj(i,j) == 1
                edge = [edge, (i,j)];
            end
        end
    end
    D = zeros(size(X,1),1);
    for i = 1:size(X,1)
        
end
function d = point_to_segment(pt, v1, v2)
    a = v1 - v2;
    b = pt - v2;
    c = pt - v1;
    if sqrt(norm(a)^2 + norm(b)^2) > norm(c)
        d = b;
    elseif sqrt(norm(a)^2 + norm(c)^2) > norm(b)
        d = c;
    else
        d = norm(cross(a,b)) / norm(a);
    end 
end

                