% network example n from slav runex1(n);
% Change theta in line 688 to adjust convolution term parameter
% Change alpha_2 in line 687 to adjust cost of travelling along network
function test_integrated(test_case, outer_iter_num, inner_iter_num)
    if test_case == 0
        Optimize_network(outer_iter_num, inner_iter_num);
    elseif test_case == 1
        Optimize_network1(outer_iter_num, inner_iter_num);
    elseif test_case == 2
        Optimize_network2(outer_iter_num, inner_iter_num);
    elseif test_case == 3
        Optimize_network3(outer_iter_num, inner_iter_num);
    elseif test_case == 4
        Optimize_network4(outer_iter_num, inner_iter_num);
    elseif test_case == 5
        Optimize_network5(outer_iter_num, inner_iter_num);
    elseif test_case == 6
        Optimize_network6(outer_iter_num, inner_iter_num);
    elseif test_case == 7
        Optimize_network7(outer_iter_num, inner_iter_num);
    elseif test_case == 8
        Optimize_network8(outer_iter_num, inner_iter_num);
    elseif test_case == 9
        Optimize_network9(outer_iter_num, inner_iter_num);
    elseif test_case == 10
        Optimize_network10(outer_iter_num, inner_iter_num);
    elseif test_case == 11
        Optimize_network11(outer_iter_num, inner_iter_num);
    elseif test_case == 12
        Optimize_network12(outer_iter_num, inner_iter_num);
    end
    
end
function Optimize_network8(outer_iter, inner_iter)
    del_r = .1;
    n_r = 12;
    x = [0,0];
    for i=1:n_r
        r_i = i*del_r;
        n_ri = round(2*pi*i);
        thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
        x = [x;r_i*cos(thetas),r_i*sin(thetas)];
    end;
    n = length(x(:,1))
    mass = 1/320*ones(1,n);
    nb = 8; rng(2);
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
    cut_indices0 = [];
    vert_indices = [];
    vert_neighs = [];
    arcs = [];
    
    lambda1 = .02;
    rho0 = .5;
    tol = 0.5*10^-4;
    alpha = .1;
    max_m = 150;
    max_avg_turn = 15;
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = 2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
            max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

function Optimize_network(outer_iter, inner_iter)
    x = [0,0;0,1;1,1;1,0];
    mass = 1/4*ones(1,4);
    [ y0,net_edges0 ] = rndInitialize(x,20);

    cut_indices0 = [];
    vert_indices = [];
    vert_neighs = [];
    arcs = [];

    lambda1 = .05;
    rho0 = .01;
    tol = 10^-4;
    alpha = .5;
    max_m = 35;
    max_avg_turn = [];
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = .2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
            max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

function Optimize_network1(outer_iter, inner_iter)
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
    y0(:,2) = zeros(m/3,1);
    y0 = [y0;zeros(m/3,1),(1/(m/3):1/(m/3):1)';zeros(m/3,1),-(1/(m/3):1/(m/3):1)'];
    cut_indices0 = [];
    vert_indices = [1,m/3,2*m/3,m];
    vert_neighs = {2,[m/3-1,m/3+1,2*m/3+1],2*m/3-1,m-1};
    arcs = {1:m/3,m/3:2*m/3,[m/3,2*m/3+1:m]};
    net_edges0 = [];

    lambda1 = .02;
    rho0 = .1;
    tol = 10^-4;
    alpha = .5;
    max_m = 30;
    max_avg_turn = [];
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = .2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
            max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end
    
function Optimize_network2(outer_iter, inner_iter)
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
    cut_indices0 = [];
    vert_indices = [1,m/3,2*m/3,m];
    vert_neighs = {2,[m/3-1,m/3+1,2*m/3+1],2*m/3-1,m-1};
    arcs = {1:m/3,m/3:2*m/3,[m/3,2*m/3+1:m]};
    net_edges0 = [];

    lambda1 = .02;
    rho0 = .2;
    tol = 10^4;
    alpha = .5;
    max_m = [];
    max_avg_turn = [];
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = .2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
        max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

function Optimize_network3(outer_iter, inner_iter)
    rng(2);
    n = 200;
    v = (0:2*pi/(n-1):2*pi)';
    x = [cos(v),sin(v)] + .1*randn(n,2);
    mass = 1/n*ones(1,n);

    m = 2*10;
    y0 = [[(0:1/(m/2-1):1)',.2*ones(m/2,1)];[flipud((0:1/(m/2-1):1)'),.2*ones(m/2,1)]];
    cut_indices0 = [];
    vert_indices = [1,m];
    vert_neighs = {[2,m],[m-1,1]};
    arcs = {1:m,[m,1]};
    net_edges0 = [];
    lambda1 = .02;
    rho0 = .1;
    tol = 10^-4;
    alpha = .5;
    max_m = 25;
    max_avg_turn = [];
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = .2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
        max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end    

function Optimize_network4(outer_iter, inner_iter)
    rng(1);
    n = 100;
    v = (0:2*pi/(n-1):2*pi)';
    x1 = [cos(v),sin(v)] + .1*randn(n,2);
    mass1 = .9/n*ones(1,n);

    x2 = [(2:3/10:5)',zeros(11,1)];
    mass2 = .1/11*ones(1,11);
    x = [x1;x2];
    mass = [mass1,mass2];

    rng(1);
    [y0,net_edges0] = rndInitialize(x,10);

    cut_indices0 = [];
    vert_indices = [];
    vert_neighs = [];
    arcs = [];

    lambda1 = .05;
    rho0 = .1;
    tol = 10^-4;
    alpha = .2;
    max_m = 25;
    max_avg_turn = [];
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = 10;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
        max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

function Optimize_network5(outer_iter, inner_iter)
    [x1,x2] = meshgrid(1:1:15,1:1:15);
    x1 = reshape(x1,length(x1)^2,1);
    x2 = reshape(x2,length(x2)^2,1);
    x = [x1,x2];
    n = length(x1);
    mass = 1/n*ones(1,n);

    rng(1);
    [ y0,net_edges0 ] = rndInitialize(x,15);

    cut_indices0 = [];
    vert_indices = [];
    vert_neighs = [];
    arcs = [];

    lambda1 = .02;
    rho0 = .01;
    tol = 10^-4;
    alpha = .5;
    max_m = 20;
    max_avg_turn = [];
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = 2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
        max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

function Optimize_network6(outer_iter, inner_iter)
    d = 3;
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
    rng(1);
    [ y0,net_edges0 ] = rndInitialize(x,30);

    cut_indices0 = [];
    vert_indices = [];
    vert_neighs = [];
    arcs = [];

    lambda1 = .02;
    rho0 = .01;
    tol = 10^-4;
    alpha = .2;
    max_m = 60;
    max_avg_turn = [];
    normalize_data = 1;
    pause_bool = 0;
    delta = .1;
    max_e_leng = .2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
        max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

function Optimize_network7(outer_iter, inner_iter)
    %del_r = r*pi/(2*(n-1))*(-1 + sqrt(1+4*(n-1)/pi)); %distance between rings, and between consecutive points on them (wrt arclength)
    %n_r = round(r/del_r);
    %del_r = r/n_r;
    del_r = .1;
    n_r = 12;
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
    mass = 1/320*ones(1,n);

    %[ y0,net_edges0 ] = rndInitialize(x,30);
    %Circle
    %         m0 = 30;
    %         thetas = (2*pi/m0:2*pi/m0:2*pi)';
    %         r0 = .7; %inner radius
    %         y0 = r0*[cos(thetas),sin(thetas)];
    %         net_edges0 = [(1:m0)',[m0,1:m0-1]'];
    %         %net_edges0 = [(2:m0)',[1:m0-1]'];
    %         %Lines
    %         nb = 4;
    %         b_len = .4;
    %         thetas = 2*pi/nb:2*pi/nb:2*pi;
    %         m = m0;
    %         for i=1:nb
    %             y0 = [y0;(r0-b_len)*cos(thetas(i)),(r0-b_len)*sin(thetas(i))];
    %             y0 = [y0;(r0+b_len)*cos(thetas(i)),(r0+b_len)*sin(thetas(i))];
    %             %y0 = [y0;r0*cos(thetas(i)),r0*sin(thetas(i))];
    %             net_edges0 = [net_edges0;round(m0*i/nb),m+1];
    %             net_edges0 = [net_edges0;round(m0*i/nb),m+2;];
    %             %m = m+ 1;
    %             m = m+2;
    %         end
    % Just Branches
    nb = 12; rng(3);
    b_len = .4;
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
        %             net_edges0 = [net_edges0;m+1,m+2];
        %             net_edges0 = [net_edges0;m+2,m+3];
        %             net_edges0 = [net_edges0;m+2,m+5];
        %             net_edges0 = [net_edges0;m+3,m+4];
        %             net_edges0 = [net_edges0;m+5,m+6];
        %                 if i<nb
        %               %     net_edges0 = [net_edges0;m+5,m+9];
        %                 end;
        m = m+6;
    end
    net_edges0 = initializeMST(y0)   ;

    %    net_edges0 = [net_edges0;3,m-1];
    %       y0 = [y0; -0.3,0; 0.3,0];
    %       net_edges0 = [net_edges0;5,9];
    %       net_edges0 = [net_edges0;21,25];
    %     net_edges0 = [net_edges0;1,17];
    %     net_edges0 = [net_edges0;12,15];
    %    net_edges0 = [net_edges0;32,35];
    %         net_edges0 = [net_edges0;8,11];
    %          net_edges0 = [net_edges0;12,15];

    % Grid


    cut_indices0 = [];
    vert_indices = [];
    vert_neighs = [];
    arcs = [];

    lambda1 = .02;
    rho0 = .5;
    tol = 0.5*10^-4;
    alpha = .1;
    max_m = 150;
    max_avg_turn = 15;
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = 2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
        max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

function Optimize_network9(outer_iter, inner_iter)
    del_r = .1;
    n_r = 7;
    x = [0,0];
    for i=1:n_r
        r_i = i*del_r;
        n_ri = round(2*pi*i);
        thetas = (2*pi/n_ri:2*pi/n_ri:2*pi)'-(pi/n_ri);
        x = [x;r_i*cos(thetas),r_i*sin(thetas)];
    end
    n = length(x(:,1))
    mass = 1/320*ones(1,n);

    nb = 3; rng(2);
    b_len = .15;
    thetas = (2*pi/nb:2*pi/nb:2*pi) + pi/nb + pi/4/nb*rand(1,nb);
    r0 = 0.3;
    y0 = [];
    m = 0;
    net_edges0 = [];
    for i=1:nb
        r0 = r0*(0.85+0.3*rand);
        y0 = [y0;(r0-b_len)*cos(thetas(i)),(r0-b_len)*sin(thetas(i))];
        y0 = [y0;r0*cos(thetas(i)),r0*sin(thetas(i))];
        y0 = [y0;(r0+ b_len)*cos(thetas(i)-pi/nb/3),(r0+ b_len)*sin(thetas(i)-pi/nb/3)]; %3
        y0 = [y0;(r0+ 2.5* b_len)*cos(thetas(i)-pi/nb/6),(r0+ 2.5* b_len)*sin(thetas(i)-pi/nb/6)]; %4
        y0 = [y0;(r0+ b_len)*cos(thetas(i)+pi/nb/3),(r0+ b_len)*sin(thetas(i)+pi/nb/3)]; % 5
        y0 = [y0;(r0+ 2.5*b_len)*cos(thetas(i)+pi/nb/6),(r0+2.5* b_len)*sin(thetas(i)+pi/nb/6)]; %6
        net_edges0 = [net_edges0;m+1,m+2];
        net_edges0 = [net_edges0;m+2,m+3];
        net_edges0 = [net_edges0;m+2,m+5];
        net_edges0 = [net_edges0;m+3,m+4];
        net_edges0 = [net_edges0;m+5,m+6];
        m = m+6;
    end
    cut_indices0 = [];
    vert_indices = [];
    vert_neighs = [];
    arcs = [];

    lambda1 = .002;
    rho0 = .5;
    tol = 0.5*10^-4;
    alpha = .1;
    max_m = 100;
    max_avg_turn = 15;
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = 2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
        max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

function Optimize_network10(outer_iter, inner_iter)
    del_r = .1;
    n_r = 17;
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
    mass = 1/320*ones(1,n);

    nb = 9; rng(1);
    b_len = .3;
    thetas = (2*pi/nb:2*pi/nb:2*pi) + pi/nb + pi/4/nb*rand(1,nb);
    r0 = 0.36;
    y0 = [0,0];
    m = 0;
    net_edges0 = [];
    for i=1:nb
        j=0;
        y0 = [y0;(r0+j*b_len)*cos(thetas(i)),(r0+j*b_len)*sin(thetas(i))];
    end
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
    end
    l=1.1*b_len;
    m=length(y0(:,1))
    for i=1:nb
        net_edges0 = [net_edges0; 1,i+1];
    end
    for i=1:m-1
        for j=i+1:m
            if norm(y0(i,:) - y0(j,:))<l
                net_edges0 = [net_edges0; i,j];
            end
        end
    end
    cut_indices0 = [];
    vert_indices = [];
    vert_neighs = [];
    arcs = [];

    lambda1 = .02;
    rho0 = .5;
    tol = 0.5*10^-4;
    alpha = .1;
    max_m = 150;
    max_avg_turn = 15;
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = 2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
        max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

function Optimize_network11(outer_iter, inner_iter)
    rng(1);
    i = (0:0.02:6.28)';
    n = length(i);
    x1 = [cos(i)+0.1*rand(n,1), sin(i)+0.1*rand(n,1)];
    x2 = x1 + [1,0];
    x=[x1;x2];
    n=2*n;

    j = (0:0.4:6)';
    m = length(j);
    l=0.5;
    y0 = [cos(j)+0.1*rand(m,1), sin(j)+0.1*rand(m,1)];
    net_edges0 = [];
    for i=1:m-1;
        for j=i+1:m
            if norm(y0(i,:) - y0(j,:))<l
                net_edges0 = [net_edges0; i,j];
            end
        end
    end
    y1 = y0+[1,0];
    y0=[y0;y1];
    m=2*m;
    % mass = ones(1,n)./n;
    net_edges1=net_edges0+[m,m];
    net_edges0=[net_edges0; net_edges1];

    mass = 1/240*ones(1,n);
    net_edges0 = [];
    for i=1:m-1;
        for j=i+1:m
            if norm(y0(i,:) - y0(j,:))<l
                net_edges0 = [net_edges0; i,j];
            end
        end
    end
    cut_indices0 = [];
    vert_indices = [];
    vert_neighs = [];
    arcs = [];

    lambda1 = .1;
    rho0 = .5;
    tol = 0.5*10^-4;
    alpha = .5;
    max_m = 120;
    max_avg_turn = 15;
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = 2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
        max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

function Optimize_network12(outer_iter, inner_iter)
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
    mass = 1/320*ones(1,n);
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
    end
    m=length(y0(:,1));
    for i=1:m-1
        for j=i+1:m
            if norm(y0(i,:) - y0(j,:))<l
                net_edges0 = [net_edges0; i,j];
            end
        end
    end
    cut_indices0 = [];
    vert_indices = [];
    vert_neighs = [];
    arcs = [];
    
    lambda1 = .02;
    rho0 = .5;
    tol = 0.5*10^-4;
    alpha = .1;
    max_m = 150;
    max_avg_turn = 15;
    normalize_data = 0;
    pause_bool = 0;
    delta = 0;
    max_e_leng = 2;
    maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
        max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n);
end

    
function maximize_subroutine(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
            max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng,outer_iter,inner_iter,n)
    for j = 1:outer_iter
        fprintf('\n outer iter = %d', j);
        tic
        [y,x,net_edges,~,~,~,~,~,~,~] = mppntest1(y0,cut_indices0,net_edges0,vert_indices,vert_neighs,arcs,x,mass,lambda1,alpha,tol,rho0,...
            max_m,max_avg_turn,normalize_data,pause_bool,delta,max_e_leng);
        toc
        p = 3;
        h = 0.005;
        alpha_2 = 0.2;
        theta = 0.5;
        sigma = 0.05;
        M = ones(n,1)/n;
        Adj = zeros(size(y,1));
        for i = 1:size(net_edges, 1)
            Adj(net_edges(i,1), net_edges(i,2)) = 1;
        end
        Adj = or(Adj, Adj');
        color_mat = rand(n,3);
        x = euler_wrap(x, y, Adj, M, inner_iter, p, h, alpha_2, theta, sigma, color_mat);
    end
end