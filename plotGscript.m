figure(1);
hold on;
for i=1:177;
    plotGeodesic(171,i,next,y,x,net_edges,edges,edge_weights,lambda,alpha);
end;
hold off;