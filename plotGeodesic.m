function [ p ] = plotGeodesic(i0,j,next,y,x,net_edges,edges,weights,lambda,alpha)
%plots the geodesic from i0 to j
% if j is empty, then the function plots the geodecis
% from i0 to all points j.
% Furthermore, i0 and j can be passed either as indices
% of the data x, or points in R^d -- in which case i0 and j
% will be se to the indices of the closest data. 
m = length(y(:,1));
V = [y;x];
N = length(V(:,1));

if length(i0)>1
    [~,i0] = min(sum((i0-x).^2,2));
end
i0 = m+i0;
rgb_c = [0,.8,0]; ls = '-';
figure;
%plotNet(x,y,net_edges,rgb_c,ls)
plotGraph(x,y,[],edges,weights,lambda,alpha,[],[],rgb_c,ls)
hold on;
p = zeros(1,length(next(:,1)));
k = 1;
i = i0;
p(k) = i;
if ~isempty(j)
    if length(j)>1
        [~,j] = min(sum((j-x).^2,2));
    end
    j = m+j;
    while next(i,j)>0
        i = next(i,j);
        vs = [V(p(k),:);V(i,:)];
        plot(vs(:,1),vs(:,2),'Color',[0,0,1],'LineWidth',.5);
        k = k+1;
        p(k) = i;
    end
    p(k:end) = [];
else
    for j = m+1:N
        k = 1;
        i = i0;
        p(k) = i;
        while next(i,j)>0
            i = next(i,j);
            vs = [V(p(k),:);V(i,:)];
            plot(vs(:,1),vs(:,2),'Color',[0,0,1],'LineWidth',.5);
            k = k+1;
            p(k) = i;
        end
        p(k:end) = [];
    end
end
scatter(V(i0,1),V(i0,2),200,[1,0,0],'.');
end
