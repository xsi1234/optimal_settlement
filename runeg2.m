function [y,net_edges,edges,edge_weights,iter,energy,next,dists,lambda,alpha] = runeg2( eg )
%RUNEG2 contains different examples of optimal networks for universal
%transport
%INPUTS:
% eg - an integer between 0 and 12 inclusive, and 51. 
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
cut_indices0 = [];
vert_indices = [];
vert_neighs = [];
arcs = [];
normalize_data = 0;
plot_bool = 1;
delta = 0;

rho0 = .2;
net_edges0 = [];  %just to initialize


switch eg
    case 0 % Example 0: Simple example with four points on corners of square
        x = [0,0;0,1;1,1;1,0];
        mass = 1/4*ones(1,4);
        
        %y0 = [0,0;.5,.5;1,1];
        %net_edges0 = [1,2;2,3];
        [ y0,net_edges0 ] = rndInitialize(x,20);
        
        lambda = .05;
        
        tol = 10^-4;
        %tol = 10^-6;
        rho0 = 1;
        alpha = .5;
        max_m = 35;
        max_avg_turn = [];
        max_e_leng = .2;
        
        
        
    case 1
        % Example 1: a horizonal 'Y'
        rng(1);
        n = 3*10;
        m = 24;
        
        x1(:,1) = -sqrt(2):sqrt(2)/(n/3-1):0;
        x1(:,2) = zeros(n/3,1);
        x2(:,1) = 0:sqrt(2)/(n/3-1):sqrt(2);
        x2(:,2) = x2(:,1);
        x3(:,1) = x2(:,1);
        x3(:,2) = -x2(:,1);
        x = [x1;x2;x3] + .1*randn(n,2);
        mass = 1/n*ones(1,n);
        
        rng(13);
        y0(:,1) = -sqrt(2):sqrt(2)/(m/3-1):0;
        %y0(:,2) = .01*rand(m,1);
        y0(:,2) = zeros(m/3,1);
        y0 = [y0;zeros(m/3,1),(1/(m/3):1/(m/3):1)';zeros(m/3,1),-(1/(m/3):1/(m/3):1)'];
        cut_indices0 = [];
        vert_indices = [1,m/3,2*m/3,m];
        vert_neighs = {2,[m/3-1,m/3+1,2*m/3+1],2*m/3-1,m-1};
        arcs = {1:m/3,m/3:2*m/3,[m/3,2*m/3+1:m]};
        
        
        lambda = .05;
        rho0=0.1;
        tol = 10^-4;
        alpha = .5;
        max_m = 30;
        max_avg_turn = [];
        max_e_leng = .2;
        
        
    case 2
        % Example 10: a horizonal 'Y'
        rng(2);
        n = 3*50;
        m = 120;
        
        x1(:,1) = -sqrt(2):sqrt(2)/(n/3-1):0;
        x1(:,2) = zeros(n/3,1);
        x2(:,1) = 0:sqrt(2)/(n/3-1):sqrt(2);
        x2(:,2) = x2(:,1);
        x3(:,1) = x2(:,1);
        x3(:,2) = -x2(:,1);
        x = [x1;x2;x3] + .1*randn(n,2);
        mass = 1/n*ones(1,n);
        
        rng(13);
        y0(:,1) = -sqrt(2):sqrt(2)/(m/3-1):0;
        y0(:,2) = zeros(m/3,1);
        y0 = [y0;zeros(m/3,1),(1/(m/3):1/(m/3):1)';zeros(m/3,1),-(1/(m/3):1/(m/3):1)'];
        vert_indices = [1,m/3,2*m/3,m];
        vert_neighs = {2,[m/3-1,m/3+1,2*m/3+1],2*m/3-1,m-1};
        arcs = {1:m/3,m/3:2*m/3,[m/3,2*m/3+1:m]};
        net_edges0 = [];
        
        lambda = .02;
        tol = 10^4;
        alpha = .5;
        max_m = [];
        max_avg_turn = [];
        max_e_leng = .2;
        
        
    case 3 % Example 3: a circle
        rng(2);
        n = 200;
        v = (0:2*pi/(n-1):2*pi)';
        x = [cos(v),sin(v)] + .1*randn(n,2);
        mass = 1/n*ones(1,n);
        
        m = 2*10;
        y0 = [[(0:1/(m/2-1):1)',.2*ones(m/2,1)];[flipud((0:1/(m/2-1):1)'),.2*ones(m/2,1)]];
        %[ y0,net_edges0 ] = rndInitialize(x,15);
        
        cut_indices0 = [];
        %vert_indices = [];
        %vert_neighs = [];
        %arcs = [];
        vert_indices = [1,m];
        vert_neighs = {[2,m],[m-1,1]};
        arcs = {1:m,[m,1]};
        net_edges0 = [];
        %vert_indices = [1,m/2,m/2+1,m];
        %vert_neighs = {[2,m],[m/2-1,m/2+1],[m/2,m/2+2],[m-1,1]};
        
        lambda = .02;
        tol = 10^-4;
        alpha = .5;
        max_m = 25;
        max_avg_turn = [];
        max_e_leng = .2;
        
    case 4 % Example 4: Circle with line segment
        rng(3);
        n = 200;
        v = (0:2*pi/(n-1):2*pi)';
        x1 = [cos(v),sin(v)] + .1*randn(n,2);
        mass1 = .9/n*ones(1,n);
        
        x2 = [(2:3/10:5)',zeros(11,1)];
        mass2 = .1/11*ones(1,11);
        x = [x1;x2];
        mass = [mass1,mass2];
        
        rng(1);
        [y0,net_edges0] = rndInitialize(x,10);
        
        
        lambda = .02;
        tol = 10^-4;
        rho0=0.4;
        alpha = .2;
        max_m = 35;
        max_avg_turn = [];
        max_e_leng = 10;
        
    case 5 %Example 5: points uniform in box
        rng(1);
        [x1,x2] = meshgrid(1:1:15,1:1:15);
        x1 = reshape(x1,length(x1)^2,1);
        x2 = reshape(x2,length(x2)^2,1);
        x = [x1,x2];
        n = length(x1);
        mass = 1/n*ones(1,n);
        m=40;
        %y0=15*rand(m,2);
        %net_edges0 = initializeMST(y0);
        [ y0,net_edges0 ] = rndInitialize(x,15);
        %lambda1 = .02;
        lambda = .07;
        tol = 10^-4;
        %alpha = .5;
        alpha = .0001;
        max_m = 120;
        max_avg_turn = [];
        max_e_leng = 2;
        
        
    case 6 %Example 6: Grid
        d = 2;
        N = 2;
        L = N+1;
        n = L*20; %points per line
        eps = .05;
        n1 = 2*n*N; %number of points from signal + noise
        n2 = 0; %number of background clutter points
        x = zeros(n1+n2,d);
        mass = 1/(n1+n2)*ones(1,n1+n2);
        rng(1);
        
        for i=1:N
            % vertical data
            x1(:,1) = i*ones(n,1) + eps*randn(n,1);
            x1(:,2) = L/n:L/n:L;
            x1(:,3:d) = eps*randn(n,d-2);
            % horizontal data
            x2(:,1) = L/n:L/n:L;
            x2(:,2) = i*ones(n,1) + eps*randn(n,1);
            x2(:,3:d) = eps*randn(n,d-2);
            x(2*(i-1)*n+1:2*i*n,:) = [x1;x2];
        end
        %x3 = [rand(n2,2)*L,zeros(n2,d-2)];
        %x3 = [rand(n2,2)*L,rand(n2,d-2)*L/2-L/4];
        %x(n1+1:n1+n2,:) = x3;
        %rng(3);
        [ y0,net_edges0 ] = rndInitialize(x,40);
%         net_edges0 = initializeMST(y0);
%         m=length(y0(:,1));
%         l=0.05;
%         for i=1:m-1;
%             for j=i+1:m
%                 if norm(y0(i,:) - y0(j,:))<l
%                     net_edges0 = [net_edges0; i,j];
%                 end
%             end
%         end;
        lambda = .02;
        alpha = .1;
        
        delta = 1; %max distance to consider adding edge
        
        rho0 = .5;
        tol = .3*10^-4;
        max_m = 100;
        max_avg_turn = [];
        normalize_data = 1;
        
        max_e_leng = .2;
        
        
        
    case 7 %Example 7:
        
        del_r = .1;
        n_r = 6;
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
        
        % Just Branches
        nb = 1; rng(2);
        b_len = .12;
        thetas = (2*pi/nb:2*pi/nb:2*pi) + pi/nb + pi/4/nb*rand(1,nb);
        r0 = 0.05;
        y0 = [];
        m = 0;
        net_edges0 = [];
        for i=1:nb
            r0 = r0*(0.85+0.3*rand);
            %  if mod(i,2)==0
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
        
        
        lambda = .02*(n/64); %n=64 <- n_r = 4
        tol = 0.5*10^-4;
        alpha = .1;
        max_m = floor(n/3);
        max_avg_turn = 15;
        max_e_leng = 2;
        
        
        
        
    case 8 %Example 8:
        
        %del_r = r*pi/(2*(n-1))*(-1 + sqrt(1+4*(n-1)/pi)); %distance between rings, and between consecutive points on them (wrt arclength)
        %n_r = round(r/del_r);
        %del_r = r/n_r;
        del_r = .1;
        n_r = 10;
        x = [0,0];
        for i=1:n_r
            r_i = i*del_r;
            n_ri = round(2*pi*i); %round(2*pi*r_i/del_r);
            %r_i = n_ri*del_r/(2*pi);
            thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
            x = [x;r_i*cos(thetas),r_i*sin(thetas)];
        end
        
        %x = r_sample.*[cos(t_sample),sin(t_sample)];
        n = length(x(:,1))
        %mass = 1/n*ones(1,n);
        mass = 1/200*ones(1,n);
        
        % Just Branches
        nb = 12; rng(1);
        b_len = .3;
        thetas = (2*pi/nb:2*pi/nb:2*pi) + pi/nb + pi/4/nb*rand(1,nb);
        r0 = 0.3;
        y0 = [];
        m = 0;
        net_edges0 = [];
        for i=1:nb
            r0 = r0*(0.2+0.8*rand);
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
        %net_edges0 = initializeMST(y0)   ;
        
        %    net_edges0 = [net_edges0;3,m-1];
        %       y0 = [y0; -0.3,0; 0.3,0];
        %       net_edges0 = [net_edges0;5,9];
        %       net_edges0 = [net_edges0;21,25];
        %     net_edges0 = [net_edges0;1,17];
        %     net_edges0 = [net_edges0;12,15];
        %    net_edges0 = [net_edges0;32,35];
        %         net_edges0 = [net_edges0;8,11];
        %          net_edges0 = [net_edges0;12,15];
        
        lambda = .02*(n/64); %n=64 <- n_r = 4
        tol = 0.5*10^-4;
        alpha = .1;
        rho0=0.3;
        %max_m = 150;
        max_m = floor(n/3);
        max_avg_turn = 15;
        max_e_leng = 2;
        
        
        
        
    case 9 % branching from middle, best for n=491, n=963
        
        del_r = .1;
        %  n_r =12; nb = 8;
        n_r = 17; nb = 10;
        rng(2);
        x = [0,0];
        for i=1:n_r
            r_i = i*del_r;
            n_ri = round(2*pi*i);
            thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
            x = [x;r_i*cos(thetas),r_i*sin(thetas)];
        end;
        n = length(x(:,1))
        mass = 1/(40*n_r)*ones(1,n);
        
        
        % Just Branches
        
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
        
        %
        %
        %
        %       %  if mod(i,2)==0
        %             y0 = [y0;(r0-b_len)*cos(thetas(i)),(r0-b_len)*sin(thetas(i))];
        %        % else
        %        %         y0 = [y0;(r0-b_len/4)*cos(thetas(i)),(r0-b_len/4)*sin(thetas(i))];
        %        % end;
        %             y0 = [y0;r0*cos(thetas(i)),r0*sin(thetas(i))];
        %             y0 = [y0;(r0+ b_len)*cos(thetas(i)-pi/nb/3),(r0+ b_len)*sin(thetas(i)-pi/nb/3)]; %3
        %             y0 = [y0;(r0+ 2.5* b_len)*cos(thetas(i)-pi/nb/6),(r0+ 2.5* b_len)*sin(thetas(i)-pi/nb/6)]; %4
        %             y0 = [y0;(r0+ b_len)*cos(thetas(i)+pi/nb/3),(r0+ b_len)*sin(thetas(i)+pi/nb/3)]; % 5
        %             y0 = [y0;(r0+ 2.5*b_len)*cos(thetas(i)+pi/nb/6),(r0+2.5* b_len)*sin(thetas(i)+pi/nb/6)]; %6
        %             y0 = [y0;(r0+ 2*b_len)*cos(thetas(i)),(r0+ 2*b_len)*sin(thetas(i))]; %7
        %         %  if mod(i,2)==0  net_edges0 = [net_edges0;m+1,m+2]; end;
        %             net_edges0 = [net_edges0;m+1,m+2]
        %             net_edges0 = [net_edges0;m+2,m+3];
        %             net_edges0 = [net_edges0;m+2,m+5];
        %             net_edges0 = [net_edges0;m+3,m+4];
        %             net_edges0 = [net_edges0;m+5,m+6];
        %             net_edges0 = [net_edges0;m+3,m+7];
        %                 if i<nb
        %                   net_edges0 = [net_edges0;m+5,m+9];
        %                 end;
        %             m = m+7;
        %         end
        %         net_edges0 = [net_edges0;3,m-1];
        %  %       y0 = [y0; -0.3,0; 0.3,0];
        %  %       net_edges0 = [net_edges0;5,9];
        %  %       net_edges0 = [net_edges0;21,25];
        %    %     net_edges0 = [net_edges0;1,17];
        %   %     net_edges0 = [net_edges0;12,15];
        %           %    net_edges0 = [net_edges0;32,35];
        % %         net_edges0 = [net_edges0;8,11];
        % %          net_edges0 = [net_edges0;12,15];
        
        lambda = .02;
        tol = 0.5*10^-4;
        alpha = .1;
        max_m = 150;
        max_avg_turn = 15;
        max_e_leng = 2;
        delta = .01;
        
        
    case 10 %RGG:
        
        %del_r = r*pi/(2*(n-1))*(-1 + sqrt(1+4*(n-1)/pi)); %distance between rings, and between consecutive points on them (wrt arclength)
        %n_r = round(r/del_r);
        %del_r = r/n_r;
        del_r = .1;
        n_r = 20;
        x = [0,0];
        for i=1:n_r
            r_i = i*del_r;
            n_ri = round(2*pi*i); %round(2*pi*r_i/del_r);
            %r_i = n_ri*del_r/(2*pi);
            thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
            x = [x;r_i*cos(thetas),r_i*sin(thetas)];
        end
        
        %x = r_sample.*[cos(t_sample),sin(t_sample)];
        n = length(x(:,1))
        %mass = 1/n*ones(1,n);
        mass = 1/(40*n_r)*ones(1,n);
        
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
        for i=1:nb;
            net_edges0 = [net_edges0; 1,i+1];
        end;
        for i=1:m-1;
            for j=i+1:m
                if norm(y0(i,:) - y0(j,:))<l
                    net_edges0 = [net_edges0; i,j];
                end
            end
        end;
        lambda = .02;
        tol = 0.5*10^-4;
        alpha = .1;
        max_m = 150;
        max_avg_turn = 15;
        max_e_leng = 2;
        
    case 11 %rand points
        
        %del_r = r*pi/(2*(n-1))*(-1 + sqrt(1+4*(n-1)/pi)); %distance between rings, and between consecutive points on them (wrt arclength)
        %n_r = round(r/del_r);
        %del_r = r/n_r;
        del_r = .1;
        n_r = 24;
        x = [0,0];
        for i=1:n_r
            r_i = i*del_r;
            n_ri = round(2*pi*i); %round(2*pi*r_i/del_r);
            %r_i = n_ri*del_r/(2*pi);
            thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
            x = [x;r_i*cos(thetas),r_i*sin(thetas)];
        end
        
        %x = r_sample.*[cos(t_sample),sin(t_sample)];
        n = length(x(:,1))
        %mass = 1/n*ones(1,n);
        mass = 1/(40*n_r)*ones(1,n);
        l=0.55;
        y0=[0,0];
        net_edges0 = [];
        rng(1);
        del_y=0.35;
        del_r = 0.3;
        r_i=0;
        for i=1:3
            r_i = r_i+del_r;
            n_ri = round(2*pi*r_i/del_y);
            %r_i = n_ri*del_r/(2*pi);
            thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
            y0 = [y0;r_i*cos(thetas),r_i*sin(thetas)];
            del_y=del_y*1.2;
            del_r=del_r*1.1;
        end;
        m=length(y0(:,1));
        for i=1:m-1;
            for j=i+1:m
                if norm(y0(i,:) - y0(j,:))<l
                    net_edges0 = [net_edges0; i,j];
                end;
            end;
        end;
        
        lambda = .02;
        tol = 0.5*10^-4;
        alpha = .1;
        max_m = 150;
        max_avg_turn = 15;
        max_e_leng = 2;
        
    case 12 %rand points
        
        %del_r = r*pi/(2*(n-1))*(-1 + sqrt(1+4*(n-1)/pi)); %distance between rings, and between consecutive points on them (wrt arclength)
        %n_r = round(r/del_r);
        %del_r = r/n_r;
        del_r = .1;
        n_r = 24;
        x = [0,0];
        for i=1:n_r
            r_i = i*del_r;
            n_ri = round(2*pi*i); %round(2*pi*r_i/del_r);
            %r_i = n_ri*del_r/(2*pi);
            thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
            x = [x;r_i*cos(thetas),r_i*sin(thetas)];
        end
        
        %x = r_sample.*[cos(t_sample),sin(t_sample)];
        n = length(x(:,1))
        %mass = 1/n*ones(1,n);
        mass = 1/(40*n_r)*ones(1,n);
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
        
        
        lambda = .02;
        tol = 0.5*10^-4;
        alpha = .1;
        max_m = 180;
        max_avg_turn = 15;
        max_e_leng = 2;
        
        
    case 51 %circle
        rng(1);
        i = (0:0.02:6.28)';
        n = length(i);
        x1 = [cos(i)+0.1*rand(n,1), sin(i)+0.1*rand(n,1)];
        x2 = x1 + [1,0];
        x=[x1;x2];
        n=2*n;
        mass = 1/240*ones(1,n);
        
%         j = (0:0.4:6)';
%         m = length(j);
%         l=0.5;
%         y0 = [cos(j)+0.1*rand(m,1), sin(j)+0.1*rand(m,1)];
%         net_edges0 = [];
%         for i=1:m-1
%             for j=i+1:m
%                 if norm(y0(i,:) - y0(j,:))<l
%                     net_edges0 = [net_edges0; i,j];
%                 end;
%             end;
%         end;
%         y1 = y0+[1,0];
%         y0=[y0;y1];
%         m=2*m;
%         % mass = ones(1,n)./n;
%         net_edges1=net_edges0+[m,m];
%         net_edges0=[net_edges0; net_edges1];
%         
%         net_edges0 = [];
%         for i=1:m-1
%             for j=i+1:m
%                 if norm(y0(i,:) - y0(j,:))<l
%                     net_edges0 = [net_edges0; i,j];
%                 end;
%             end;
%         end;
        
        [ y0,net_edges0 ] = rndInitialize(x,30);
        
        lambda = .4;
        tol = 0.5*10^-4;
        alpha = .001;
        max_m = 120;
        max_avg_turn = 15;
        max_e_leng = 2;
        %rho = 1;
end

if isempty(net_edges0)
    m = length(y0(:,1));
    neighbors = cell(m,1); % cell containing the indices corresponding to the neighbors of each point
    neighbors(vert_indices) = vert_neighs; %elements of neighbors need not be
    % ordered, since net_edges need not be oriented in any particular direction.
    arc_inds = cell(m,1); %entry i gives the arc indices that y_i is contained in
    y_net_edge_inds = cell(m,1); %entry i gives the net_edge indices that y_i neighbors
    y_verts = cell(m,1); %entry i gives the topological vert indices that y_i 'neighbors'
    
    for i=1:length(arcs)
        arc_inds{arcs{i}(1)} = [arc_inds{arcs{i}(1)},i];
        y_verts{arcs{i}(1)} = [y_verts{arcs{i}(1)},arcs{i}(end)];
        for j=2:length(arcs{i})-1
            neighbors{arcs{i}(j)} = [arcs{i}(j-1),arcs{i}(j+1)];
            arc_inds{arcs{i}(j)} = [arc_inds{arcs{i}(j)},i];
            y_verts{arcs{i}(j)} = [y_verts{arcs{i}(j)},setdiff([arcs{i}(1),arcs{i}(end)],arcs{i}(j))];
        end
        arc_inds{arcs{i}(end)} = [arc_inds{arcs{i}(end)},i];
        y_verts{arcs{i}(end)} = [y_verts{arcs{i}(end)},arcs{i}(1)];
    end
    
    num_net_edges = m-length(vert_indices) + length(cell2mat(vert_neighs))/2.0;
    net_edges0 = zeros(num_net_edges,2);
    k = 1;
    for i=1:m
        for j=1:length(neighbors{i})
            if neighbors{i}(j)>i
                net_edges0(k,:) = [i,neighbors{i}(j)];
                y_net_edge_inds{i} = [y_net_edge_inds{i},k];
                y_net_edge_inds{neighbors{i}(j)} = [y_net_edge_inds{neighbors{i}(j)},k];
                k = k + 1;
            end
        end
    end
end


tic;
[y,net_edges,edges,edge_weights,iter,energy,dists,next] = onut(y0,net_edges0,x,mass,lambda,alpha,...
    tol,rho0,max_m,max_avg_turn,normalize_data,plot_bool,delta,max_e_leng);
toc;
%if exist('n')
%    fname = ['results2/n',num2str(n)];
%    save(fname)
%end


end

