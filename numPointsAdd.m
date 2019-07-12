function n_add = numPointsAdd(y,max_avg_turn,max_m,cut_indices)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

m = length(y(:,1));

if ~isempty(max_avg_turn)
    v = y(2:m,:)-y(1:m-1,:);
    normv = sqrt(sum(v.^2,2));
    vdots = sum(bsxfun(@times,v(2:m-1,:),v(1:m-2,:)),2);
    ta_deg = acosd(bsxfun(@rdivide,vdots,bsxfun(@times,normv(1:m-2),normv(2:m-1))));
    comp_ends = union(cut_indices(cut_indices>1)-1,cut_indices(cut_indices<m-1));
    ta_deg(comp_ends) = 0;
    
    avg_turn_deg = sum(ta_deg)/(m - 2 - length(comp_ends));
    
    if avg_turn_deg > max_avg_turn
        n_add = max(1,floor(m * min(.4, sqrt(avg_turn_deg/max_avg_turn - 1))));
        n_add = min(max_m-m,n_add);
    else
        n_add = 0;
    end
    start_ind = [1;cut_indices-1];
    end_ind = [cut_indices;m];
    %n_add = min(n_add+sum(end_ind == start_ind + 1),floor(.4*m)); %In case there are components with length 2 (no turning angle)
    %if ~isempty(max_m)
    
    %end
else
    n_add = max(min(max_m-m,max(1,floor(.5*m))),0);
end

end

