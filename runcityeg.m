function [y,x,net_edges,edges,edge_weights,iter,energy,dists,next,lambda,alpha] = runcityeg( eg )
%RUNCITYEG examples with different initializations for distribution in a ball
%   The function runs onut.m with data uniformly distributed in a ball, and
%   initialization given by the specified example number below. The size of
%   the city is controlled by the n_r (number of rings) parameter below,
%   along with other parameters that are set globally. 
%INPUTS:
% eg - an integer between -1 and 12 inclusive, but not including 6,
%   specifying which initialization to use. 
%OUTPUTS:
% y - locations of points on final network, given as mxd vector
%   (m-number of points, d-dimension)
% x - the n data points, given as nxd vector
% net_edges - a list of edges for the initial network. Each edge contains
%   the indices of the points in y0 that form the edge. It is given as a
%   n_e x 2 vector, where n_e is the number of edges
% edges - a list of edges on the graph (XuY), whose vertices consist of all
%   data points 'x', and all network points 'y'.
% edge_weights - a list of weights corresponding to 'edges' (of the graph
%   XuY). This vector contains the total mass traveling through each edge,
%   'm_{a,b}' as defined in the paper.
% iter - the total number of (outer) iterations of the onut algorithm.
% energy - the discrete energy of the final network (y,net_edges).
% dists - the geodesic (shortest path) distances between all m+n vertices 
%   on the complete graph XuY, with weights equal to length of the edge,
%   times alpha if the edge belongs to the network.
% next - an (m+n)x(m+n) matrix of indices indicating the shortest paths
%   between all pairs of vertices on the graph XuY. the entry at next(i,j)
%   gives the point that should be visited next if going from point i to
%   point j. 
% lambda - parameter in ONST functional (value > 0)
% alpha - cost of travelling along network, value between 0 and 1.

%universal parameters
normalize_data = 0;
plot_bool = 1;
delta = .7;

n_r = 20;
del_r = .1;
x = [0,0];
for i=1:n_r
    r_i = i*del_r;
    n_ri = round(2*pi*i); %round(2*pi*r_i/del_r);
    %r_i = n_ri*del_r/(2*pi);
    thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
    x = [x;r_i*cos(thetas),r_i*sin(thetas)];
end
n = length(x(:,1))
%mass = 1/(50*n_r)*ones(1,n);
mass = 1/200*ones(1,n);

%rho0 = .2;
rho0 = n/200; %should increase as the spacing between points (on network) decreases
lambda = .02*(n/64); %n=64 <- n_r = 4
%tol = 0.5*10^-4;
tol = 10^-5;
alpha = .1;
max_m = floor(n/2);
max_avg_turn = 15;
max_e_leng = del_r;

net_edges0 = [];
y0 = [];

switch eg
    case -1 %line segment
        %y0 = [.01,0;-.01,0];
        %net_edges0 = [1,2];
        r = .2;
        len = .4;
        seg = [0,len/2;0,-len/2];
        %multiple segments rotated
        y0 = [];
        net_edges0 = [];
        k = 3;
        theta = (2*pi/k:2*pi/k:2*pi);
        for i=1:length(theta)
            rot_mat = [[cos(theta(i)),-sin(theta(i))];[sin(theta(i)),cos(theta(i))]];
            new_seg = r*[cos(theta(i)),sin(theta(i))] + (rot_mat*seg')';
            y0 = [y0;new_seg];
            net_edges0 = [net_edges0;2*i-1,2*i];
        end
        
        
    case 0
        %[ y0,net_edges0 ] = rndInitialize(x,100);
        rng(1);
        %y0 = datasample(x,200,'Replace',false);
        y0 = 2*[-.1,.1;.1,.1;0,.03;0,0;0,-.03;-.1,-.1;.1,-.1]; %best for n_r = 5;
        %y0 = 3*[-.1,.1;.1,.1;0,.03;0,0;0,-.03;-.1,-.1;.1,-.1]; 
        net_edges0 = initializeMST(y0);
        
    case 1 %Circle

        m0 = 30;
        thetas = (2*pi/m0:2*pi/m0:2*pi)';
        r0 = .6; %inner radius
        y0 = r0*[cos(thetas),sin(thetas)];
        net_edges0 = [(1:m0)',[m0,1:m0-1]'];
        
        nb = 10;
        b_len = 1; 
        r0 = .2;
        thetas = 2*pi/nb:2*pi/nb:2*pi;
        m = m0;
        for i=1:nb
            %y0 = [y0;(r0-b_len)*cos(thetas(i)),(r0-b_len)*sin(thetas(i))];
            %y0 = [y0;(r0+b_len)*cos(thetas(i)),(r0+b_len)*sin(thetas(i))];
            y0 = [y0;r0*cos(thetas(i)),r0*sin(thetas(i))];
            y0 = [y0;(r0+b_len)*cos(thetas(i)),(r0+b_len)*sin(thetas(i))];
            %y0 = [y0;r0*cos(thetas(i)),r0*sin(thetas(i))];
            net_edges0 = [net_edges0;round(m0*i/nb),m+1];
            net_edges0 = [net_edges0;round(m0*i/nb),m+2;];
            %m = m+ 1;
            m = m+2;
        end
        
    case 2 %Simple branches
        
        nb = 3;
        b_len = .1;
        thetas = (2*pi/nb:2*pi/nb:2*i) + pi/nb;
        r0 = .2; %inner radius
        y0 = [0,0];
        m = 1;
        net_edges0 = [];
        for i=1:nb
            y0 = [y0;(r0-b_len)*cos(thetas(i)),(r0-b_len)*sin(thetas(i))];
            %y0 = [y0;r0*cos(thetas(i)),r0*sin(thetas(i))];
            y0 = [y0;(r0+b_len)*cos(thetas(i)),(r0+b_len)*sin(thetas(i))];
            net_edges0 = [net_edges0;1,m+1;m+1,m+2];
            m = m+2;
        end
        
    
    case 3 % Tree
        nl = 2; %number of layers
        nb = 3; %number of branches
        r = .7; %scale factor
        b_len = .5; %initial branch length
        thetas = (2*pi/nb:2*pi/nb:2*pi);
        y0 = [0,0];
        m = 1;
        net_edges0 = [];
        end_pts = [0,0];
        end_theta = 0;
        end_inds = 1;
 
        for i=1:nl
            new_end_pts = [];
            new_end_theta = [];
            new_end_inds = [];
            for j=1:length(end_pts(:,1))
                for k=1:nb
                    new_theta = end_theta(j)+thetas(k);
                    new_point = [end_pts(j,1)+b_len*cos(new_theta),end_pts(j,2)+b_len*sin(new_theta)];
                    y0 = [y0;new_point];
                    m = m+1;
                    net_edges0 = [net_edges0;end_inds(j),m];
                    new_end_pts = [new_end_pts;new_point];
                    new_end_theta = [new_end_theta; new_theta];
                    new_end_inds = [new_end_inds;m];
                end
            end
            end_pts = new_end_pts;
            end_theta = new_end_theta;
            end_inds = new_end_inds;
            b_len = b_len*r;
        end
        
        
    case 4 %tree where all vertices are deg 3
        rng(1);
        %nl = 3; %number of layers
        nl = 3;
        nb = 4; %number of initial branches
        r = .7; %scale factor
        %r = 1.6;
        %b_len = .6; %initial branch length
        b_len = .45;
        thetas = (2*pi/nb:2*pi/nb:2*pi);
        y0 = [0,0];
        m = 1;
        net_edges0 = [];
        end_pts = [0,0];
        end_theta = 0;
        end_inds = 1;
        new_end_pts = [];
        new_end_theta = [];
        new_end_inds = [];
        for k=1:nb
            new_theta = end_theta+thetas(k);
            new_point = [end_pts(1,1)+b_len*cos(new_theta),end_pts(1,2)+b_len*sin(new_theta)];
            y0 = [y0;new_point];
            m = m+1;
            net_edges0 = [net_edges0;end_inds,m];
            new_end_pts = [new_end_pts;new_point];
            new_end_theta = [new_end_theta; new_theta];
            new_end_inds = [new_end_inds;m];
        end
        end_pts = new_end_pts;
        end_theta = new_end_theta;
        end_inds = new_end_inds;
        b_len = b_len*r;
        nb1 = 2;
        thetas = [pi/3,-pi/3];
        %thetas = [pi/2,-pi/2];
        %thetas = [pi/6,-pi/6];
        for i=1:nl
            new_end_pts = [];
            new_end_theta = [];
            new_end_inds = [];
            for j=1:length(end_pts(:,1))
                for k=1:nb1
                    new_theta = end_theta(j)+thetas(k);
                    new_point = [end_pts(j,1)+b_len*cos(new_theta),end_pts(j,2)+b_len*sin(new_theta)];
                    y0 = [y0;new_point];
                    m = m+1;
                    net_edges0 = [net_edges0;end_inds(j),m];
                    new_end_pts = [new_end_pts;new_point];
                    new_end_theta = [new_end_theta; new_theta];
                    new_end_inds = [new_end_inds;m];
                end
            end
            end_pts = new_end_pts;
            end_theta = new_end_theta;
            end_inds = new_end_inds;
            b_len = b_len*r;
        end
        
        r = .25;
        len = .2;
        seg = [0,len;0,0];
        %multiple segments rotated
        y1 = [];
        net_edges1 = [];
        k = nb;
        theta = (2*pi/k:2*pi/k:2*pi);
        for i=1:length(theta)
            rot_mat = [[cos(theta(i)),-sin(theta(i))];[sin(theta(i)),cos(theta(i))]];
            new_seg = r*[cos(theta(i)),sin(theta(i))] + (rot_mat*seg')';
            y1 = [y1;new_seg];
            net_edges1 = [net_edges1;2*i-1,2*i];
        end
        
        m = length(y0(:,1));
        y0 = 1.4*[y0;y1];
        net_edges0 = [net_edges0;net_edges1+m];
        
        
    case 5
        
        rng(1);
        nl = 1; %number of layers
        nb = 3; %number of branches
        r = 1.5; %scale factor
        b_len = .5; %initial branch length
        thetas = (2*pi/nb:2*pi/nb:2*pi);
        y0 = [y0;0,0];
        m = length(y0(:,1));
        %net_edges0 = [];
        end_pts = [0,0];
        end_theta = 0;
        end_inds = m;
        new_end_pts = [];
        new_end_theta = [];
        new_end_inds = [];
        for k=1:nb
            new_theta = end_theta+thetas(k);
            new_point = [end_pts(1,1)+b_len*cos(new_theta),end_pts(1,2)+b_len*sin(new_theta)];
            y0 = [y0;new_point];
            m = m+1;
            net_edges0 = [net_edges0;end_inds,m];
            new_end_pts = [new_end_pts;new_point];
            new_end_theta = [new_end_theta; new_theta];
            new_end_inds = [new_end_inds;m];
        end
        end_pts = new_end_pts;
        end_theta = new_end_theta;
        end_inds = new_end_inds;
        b_len = b_len*r;
        nb = 2;
        %thetas = [pi/3,-pi/3];
        thetas = [2*pi/5,-2*pi/5];
        for i=1:nl
            new_end_pts = [];
            new_end_theta = [];
            new_end_inds = [];
            for j=1:length(end_pts(:,1))
                for k=1:nb
                    new_theta = end_theta(j)+thetas(k);
                    new_point = [end_pts(j,1)+b_len*cos(new_theta),end_pts(j,2)+b_len*sin(new_theta)];
                    y0 = [y0;new_point];
                    m = m+1;
                    net_edges0 = [net_edges0;end_inds(j),m];
                    new_end_pts = [new_end_pts;new_point];
                    new_end_theta = [new_end_theta; new_theta];
                    new_end_inds = [new_end_inds;m];
                end
            end
            end_pts = new_end_pts;
            end_theta = new_end_theta;
            end_inds = new_end_inds;
            b_len = b_len*r;
        end
        %circle
        %r0 = norm(end_pts(1,:)); %radius
        r0 = max(sqrt(sum(end_pts.^2,2)));
        m0 = 30;
        thetas = (2*pi/m0:2*pi/m0:2*pi)';
        y1 = r0*[cos(thetas),sin(thetas)];
        net_edges1 = [(1:m0)',[m0,1:m0-1]'];
        net_edges0 = [net_edges0;net_edges1+length(y0(:,1))];
        y0 = [y0;y1];
        
        %second layer of branches
        deg = 2;
        nl = 2;
        nb = 12;
        b_len = .5;
        r = .4;
        %thetas = [pi/3,-pi/3];
        thetas = [pi/6,-pi/6];
        m = length(y0(:,1));
        end_theta = (2*pi/nb:2*pi/nb:2*pi);
        end_pts = r0*[cos(end_theta)',sin(end_theta)'];
        y0 = [y0;end_pts];
        end_inds = m+1:m+length(end_pts(:,1));
        m = length(y0(:,1));
        new_end_pts = [];
        %new_end_theta = [];
        new_end_inds = [];
        
        for k=1:nb
            %new_theta = end_theta;
            new_point = [end_pts(k,1)+b_len*cos(end_theta(k)),end_pts(k,2)+b_len*sin(end_theta(k))];
            y0 = [y0;new_point];
            m = m+1;
            net_edges0 = [net_edges0;end_inds(k),m];
            new_end_pts = [new_end_pts;new_point];
            %new_end_theta = [new_end_theta; new_theta];
            new_end_inds = [new_end_inds;m];
        end
        end_pts = new_end_pts;
        %end_theta = new_end_theta;
        end_inds = new_end_inds;
        
        for i=1:nl
            new_end_pts = [];
            new_end_theta = [];
            new_end_inds = [];
            for j=1:length(end_pts(:,1))
                for k=1:deg
                    new_theta = end_theta(j)+thetas(k);
                    new_point = [end_pts(j,1)+b_len*cos(new_theta),end_pts(j,2)+b_len*sin(new_theta)];
                    y0 = [y0;new_point];
                    m = m+1;
                    net_edges0 = [net_edges0;end_inds(j),m];
                    new_end_pts = [new_end_pts;new_point];
                    new_end_theta = [new_end_theta; new_theta];
                    new_end_inds = [new_end_inds;m];
                end
            end
            end_pts = new_end_pts;
            end_theta = new_end_theta;
            end_inds = new_end_inds;
            b_len = b_len*r;
        end
        y0 = .8*y0;
        
    case 7 %2 layers of branches
        
        % Just Branches
        nb = 9; rng(2);
        b_len = .08;
        thetas = (2*pi/nb:2*pi/nb:2*pi) + pi/nb + pi/4/nb*rand(1,nb);
        r0 = 0.15;
        y0 = [];
        m = 0;
        net_edges0 = [];
        for i=1:nb
            %r0 = r0*(0.85+0.3*rand);
            %r0 = .7*r0;
            %b_len = .7*b_len;
            %  if mod(i,2)==0
            y0 = [y0;(r0-b_len)*cos(thetas(i)),(r0-b_len)*sin(thetas(i))];
            % else
            % y0 = [y0;(r0-b_len/4)*cos(thetas(i)),(r0-b_len/4)*sin(thetas(i))];
            % end;
            y0 = [y0;r0*cos(thetas(i)),r0*sin(thetas(i))];
            y0 = [y0;(r0+ b_len)*cos(thetas(i)-pi/nb/3),(r0+ b_len)*sin(thetas(i)-pi/nb/3)]; %3
            y0 = [y0;(r0+ 2.5* b_len)*cos(thetas(i)-pi/nb/4),(r0+ 2.5* b_len)*sin(thetas(i)-pi/nb/4)]; %4
            y0 = [y0;(r0+ b_len)*cos(thetas(i)+pi/nb/3),(r0+ b_len)*sin(thetas(i)+pi/nb/3)]; % 5
            y0 = [y0;(r0+ 2.5*b_len)*cos(thetas(i)+pi/nb/4),(r0+2.5* b_len)*sin(thetas(i)+pi/nb/4)]; %6
            %  y0 = [y0;(r0+ 2*b_len)*cos(thetas(i)),(r0+ 2*b_len)*sin(thetas(i))]; %7
            %  if mod(i,2)==0  net_edges0 = [net_edges0;m+1,m+2]; end;
            net_edges0 = [net_edges0;m+1,m+2];
            net_edges0 = [net_edges0;m+2,m+3];
            net_edges0 = [net_edges0;m+2,m+5];
            net_edges0 = [net_edges0;m+3,m+4];
            net_edges0 = [net_edges0;m+5,m+6];
            %    net_edges0 = [net_edges0;m+3,m+7];
            %  if i<nb
            %   net_edges0 = [net_edges0;m+5,m+9];
            % end;
            m = m+6;
        end
        %      net_edges0 = [net_edges0;3,m-1];
        %       y0 = [y0; -0.3,0; 0.3,0];
        %       net_edges0 = [net_edges0;5,9];
        %       net_edges0 = [net_edges0;21,25];
        %     net_edges0 = [net_edges0;1,17];
        %     net_edges0 = [net_edges0;12,15];
        %    net_edges0 = [net_edges0;32,35];
        %         net_edges0 = [net_edges0;8,11];
        %          net_edges0 = [net_edges0;12,15];
        y0 = 4*y0;
        
    case 8
                % Just Branches
        nb = 5; rng(1);
        b_len = .3;
        thetas = (2*pi/nb:2*pi/nb:2*pi) + pi/nb + pi/4/nb*rand(1,nb);
        r0 = 0.5;
        y0 = [];
        m = 0;
        net_edges0 = [];
        for i=1:nb
            %r0 = r0*(0.2+0.8*rand);
            y0 = [y0;(r0-b_len)*cos(thetas(i)),(r0-b_len)*sin(thetas(i))];
            y0 = [y0;r0*cos(thetas(i)),r0*sin(thetas(i))];
            y0 = [y0;(r0+ 0.8*b_len)*cos(thetas(i)-pi/nb/2.5),(r0+ 0.8*b_len)*sin(thetas(i)-pi/nb/2.5)]; %3
            y0 = [y0;(r0+ 1.4* b_len)*cos(thetas(i)-pi/nb/3),(r0+ 1.4* b_len)*sin(thetas(i)-pi/nb/3)]; %4
            y0 = [y0;(r0+ 0.8*b_len)*cos(thetas(i)+pi/nb/2.5),(r0+ 0.8*b_len)*sin(thetas(i)+pi/nb/2.5)]; % 5
            y0 = [y0;(r0+ 1.4*b_len)*cos(thetas(i)+pi/nb/3),(r0+1.4* b_len)*sin(thetas(i)+pi/nb/3)]; %6
            net_edges0 = [net_edges0;m+1,m+2];
            net_edges0 = [net_edges0;m+2,m+3];
            net_edges0 = [net_edges0;m+2,m+5];
            net_edges0 = [net_edges0;m+3,m+4];
            net_edges0 = [net_edges0;m+5,m+6];
            if i<nb
                net_edges0 = [net_edges0;m+5,m+9];
            end;
            m = m+6;
        end
        
        
        
    case 9 %Ring with longe wedges

        nb = 10;
        b_len = .3;
        thetas = (2*pi/nb:2*pi/nb:2*pi) + pi/nb + pi/4/nb*rand(1,nb);
        r0 = 0.4;
        y0 = [];
        m = 0;
        net_edges0 = [];
        for i=1:nb
            r0 = r0*(0.85+0.3*rand);
            if mod(i,2)==0
                y0 = [y0;(r0-b_len)*cos(thetas(i)),(r0-b_len)*sin(thetas(i))];
            else
                y0 = [y0;(r0-b_len/4)*cos(thetas(i)),(r0-b_len/4)*sin(thetas(i))];
            end;
            y0 = [y0;r0*cos(thetas(i)),r0*sin(thetas(i))];
            y0 = [y0;(r0+ b_len)*cos(thetas(i)-pi/nb/3),(r0+ b_len)*sin(thetas(i)-pi/nb/3)]; %3
            y0 = [y0;(r0+ 2.5* b_len)*cos(thetas(i)-pi/nb/6),(r0+ 2.5* b_len)*sin(thetas(i)-pi/nb/6)]; %4
            y0 = [y0;(r0+ b_len)*cos(thetas(i)+pi/nb/3),(r0+ b_len)*sin(thetas(i)+pi/nb/3)]; % 5
            y0 = [y0;(r0+ 2.5*b_len)*cos(thetas(i)+pi/nb/6),(r0+2.5* b_len)*sin(thetas(i)+pi/nb/6)]; %6
            %  y0 = [y0;(r0+ 2*b_len)*cos(thetas(i)),(r0+ 2*b_len)*sin(thetas(i))]; %7
            if mod(i,2)==0  net_edges0 = [net_edges0;m+1,m+2]; end;
            %   net_edges0 = [net_edges0;m+1,m+2]
            net_edges0 = [net_edges0;m+2,m+3];
            net_edges0 = [net_edges0;m+2,m+5];
            net_edges0 = [net_edges0;m+3,m+4];
            net_edges0 = [net_edges0;m+5,m+6];
            %    net_edges0 = [net_edges0;m+3,m+7];
            if i<nb
                net_edges0 = [net_edges0;m+5,m+9];
            end;
            m = m+6;
        end
        net_edges0 = [net_edges0;3,m-1];
        
        
        
    case 10 %Branches with ring and triangles
        nb = 9; rng(2);
        b_len = .34;
        thetas = (2*pi/nb:2*pi/nb:2*pi) + pi/nb + pi/4/nb*rand(1,nb);
        r0 = 0.42;
        y0 = [0,0];
        m = 0;
        net_edges0 = [];
        for i=1:nb
            j=0;
            y0 = [y0;(r0+j*b_len)*cos(thetas(i)),(r0+j*b_len)*sin(thetas(i))];
        end;
        for i=1:nb
            j=0.6;
            y0 = [y0;(r0+j*b_len)*cos(thetas(i)),(r0+j*b_len)*sin(thetas(i))];
            j=2;
            y0 = [y0;(r0+j*b_len)*cos(thetas(i)),(r0+j*b_len)*sin(thetas(i))];
        end
        j=2.9;
        for i=1:nb
            y0 = [y0;(r0+j*b_len)*cos(thetas(i)+pi/nb/4),(r0+j*b_len)*sin(thetas(i)+pi/nb/4)];
            y0 = [y0;(r0+j*b_len)*cos(thetas(i)-pi/nb/4),(r0+j*b_len)*sin(thetas(i)-pi/nb/4)];
        end;
        l=1.1*b_len;
        m=length(y0(:,1))
        for i=1:nb
            net_edges0 = [net_edges0; 1,i+1];
        end;
        for i=1:m-1
            for j=i+1:m
                if norm(y0(i,:) - y0(j,:))<l
                    net_edges0 = [net_edges0; i,j];
                end
            end
        end;
        
        
    case 11 % mesh
        
        l=0.55;
        y0=[0,0];
        net_edges0 = [];
        rng(1);
        del_y=0.35;
        del_r = 0.3;
        r_i=0;
        for i=1:5
            r_i = r_i+del_r;
            n_ri = round(2*pi*r_i/del_y);
            %r_i = n_ri*del_r/(2*pi);
            thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
            y0 = [y0;r_i*cos(thetas),r_i*sin(thetas)];
            del_y=del_y*1.2;
            del_r=del_r*1.1;
        end;
        m=length(y0(:,1));
        for i=1:m-1
            for j=i+1:m
                if norm(y0(i,:) - y0(j,:))<l
                    net_edges0 = [net_edges0; i,j];
                end;
            end;
        end;
        
        
    case 12 %Rings with wedges
        
        l=0.85;
        y0=[0,0];
        net_edges0 = [];
        rng(1);
        del_y=0.6;
        del_r = 0.65;
        r_i=0;
        for i=1:2
            r_i = r_i+del_r;
            n_ri = round(2*pi*r_i/del_y);
            %r_i = n_ri*del_r/(2*pi);
            thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
            y0 = [y0;r_i*cos(thetas),r_i*sin(thetas)];
            del_y=del_y*1.2;
            del_r=del_r*1.1;
        end;
        m=length(y0(:,1));
        for i=1:m-1
            for j=i+1:m
                if norm(y0(i,:) - y0(j,:))<l
                    net_edges0 = [net_edges0; i,j];
                end;
            end;
        end;
        
        
        nb=14;
        b_len = .4;
        thetas = (2*pi/nb:2*pi/nb:2*pi) + pi/nb + pi/4/nb*rand(1,nb);
        r0 = 1.4;
        % y0 = [];
        % m = 0;
        % net_edges0 = [];
        for i=1:nb
            r0 = r0*(0.9+0.2*rand);
            % if mod(i,2)==0
            y0 = [y0;(r0-b_len)*cos(thetas(i)),(r0-b_len)*sin(thetas(i))];
            % else
            %         y0 = [y0;(r0-b_len/4)*cos(thetas(i)),(r0-b_len/4)*sin(thetas(i))];
            % end;
            y0 = [y0;r0*cos(thetas(i)),r0*sin(thetas(i))];
            y0 = [y0;(r0+ b_len)*cos(thetas(i)-pi/nb/3),(r0+ b_len)*sin(thetas(i)-pi/nb/3)]; %3
            y0 = [y0;(r0+ 2.5* b_len)*cos(thetas(i)-pi/nb/4),(r0+ 2.5* b_len)*sin(thetas(i)-pi/nb/4)]; %4
            y0 = [y0;(r0+ b_len)*cos(thetas(i)+pi/nb/3),(r0+ b_len)*sin(thetas(i)+pi/nb/3)]; % 5
            y0 = [y0;(r0+ 2.5*b_len)*cos(thetas(i)+pi/nb/4),(r0+2.5* b_len)*sin(thetas(i)+pi/nb/4)]; %6
            %  y0 = [y0;(r0+ 2*b_len)*cos(thetas(i)),(r0+ 2*b_len)*sin(thetas(i))]; %7
            % if mod(i,2)==0  net_edges0 = [net_edges0;m+1,m+2]; end;
            net_edges0 = [net_edges0;m+1,m+2];
            net_edges0 = [net_edges0;m+2,m+3];
            net_edges0 = [net_edges0;m+2,m+5];
            net_edges0 = [net_edges0;m+3,m+4];
            net_edges0 = [net_edges0;m+5,m+6];
            %    net_edges0 = [net_edges0;m+3,m+7];
            if i<nb
                net_edges0 = [net_edges0;m+5,m+9];
            end;
            m = m+6;
        end
        net_edges0 = [net_edges0;3,m-1];
        
        
end

tic;
[y,net_edges,edges,edge_weights,iter,energy,dists,next] = onut(y0,net_edges0,x,mass,lambda,alpha,...
    tol,rho0,max_m,max_avg_turn,normalize_data,plot_bool,delta,max_e_leng);
toc;

%save file and fig
m = length(y(:,1));
es = num2str(energy);
es = es(es ~= '.');
es_sig = find(es ~= '0');
es_st = es_sig(1); %index of starting digit
fname = ['results2/n',num2str(n),'e',es(es_st:es_st+3),'m',num2str(m),'ex',num2str(eg)];
while exist(fname,'file') ~= 0
    fname = [fname,'&'];
end
$save(fname);
%savefig(fname);

end

