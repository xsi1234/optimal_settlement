function [ ] = plotNet( x,y,net_edges,rgb_c,ls)
%PLOTNET Plots the data and the current network
%   The edges of the transportation network y are plotted in color 'rgb_c',
%   with linestyle 'ls'.

%plotADP(x,y,arcs,mass,[],0,rgb_c,ls); 
scatter(x(:,1),x(:,2),40,[0,0,0],'.');
hold on;
for i=1:length(net_edges(:,1))
    plot(y(net_edges(i,:),1),y(net_edges(i,:),2),'Color',rgb_c,'LineWidth',3,'LineStyle',ls);
end
scatter(y(:,1),y(:,2),4,'red');
% n_edges = edges(edges(:,1)<=m & edges(:,2)<=m & edge_weights<lambda/(1-alpha),:);
% for i=1:length(n_edges(:,1))
% plot([y(n_edges(i,1),1);y(n_edges(i,2),1)],[y(n_edges(i,1),2);y(n_edges(i,2),2)],'blue');
% end
% nd_edges = edges(edges(:,1)>m & edges(:,2)<=m,:);
% for i=1:length(nd_edges(:,1))
% plot([x(nd_edges(i,1)-m,1);y(nd_edges(i,2),1)],[x(nd_edges(i,1)-m,2);y(nd_edges(i,2),2)],'blue');
% end
% data_edges = edges(edges(:,1)>m & edges(:,2)>m,:);
% for i=1:length(data_edges(:,1))
% plot(x(data_edges(i,:)-m,1),x(data_edges(i,:)-m,2),'blue');
% end
hold off;
axis equal;


end
