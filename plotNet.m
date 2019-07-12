function [ ] = plotNet( x,y,net_edges,rgb_c,ls)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%plotADP(x,y,arcs,mass,[],0,rgb_c,ls); 
scatter(x(:,1),x(:,2),40,[0,0,0],'.');
hold on;
for i=1:length(net_edges(:,1))
    plot(y(net_edges(i,:),1),y(net_edges(i,:),2),'Color',rgb_c,'LineWidth',3,'LineStyle',ls);
end
scatter(y(:,1),y(:,2),4,'red');
hold off;
axis equal;


end
