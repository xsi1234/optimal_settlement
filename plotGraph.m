function [ ] = plotGraph( x,y,edges,weights,lambda,alpha,num,mass,rgb_c,ls)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% cmap = colormap(jet);
% clen = length(cmap(:,1));
% cmap = cmap(floor(clen/2):end,:);
% clen = length(cmap(:,1));
%plotADP(x,y,arcs,mass,[],0,rgb_c,ls);
scatter(x(:,1),x(:,2),45,[0,0,0],'.');
hold on;
v = [y;x];
max_w = max(weights);
[~,sort_ind] = sort(weights,'descend');
k = find(weights(sort_ind)>lambda/(1-alpha),1,'last');
min_w = weights(sort_ind(k));
if max_w == min_w
    for i=1:k
        ind = edges(sort_ind(i),:);
        plot(v(ind,1),v(ind,2),'Color',[1,0,0],'LineWidth',weights(sort_ind(i))*50,'LineStyle',ls);
        %plot(v(ind,1),v(ind,2),'Color',cmap(end,:),'LineWidth',weights(sort_ind(i))*50,'LineStyle',ls);
    end
else
    for i=1:k
        ind = edges(sort_ind(i),:);
        plot(v(ind,1),v(ind,2),'Color',[(max_w-weights(sort_ind(i)))/(max_w-min_w),0,0],'LineWidth',1.5+weights(sort_ind(i))*3/max_w,'LineStyle',ls);
        %plot(v(ind,1),v(ind,2),'Color',cmap(max(1,round((weights(sort_ind(i))-min_w)/(max_w-min_w)*clen)),:),'LineWidth',weights(sort_ind(i))*50,'LineStyle',ls);
    end
end
hold off;
axis equal;


end

