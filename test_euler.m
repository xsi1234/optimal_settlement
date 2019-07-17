%alpha: proportion of cost of travelling along the network to travelling directly
%theta: repulsion parameter
%sigma: variance in K
%h: step size for forward Euler Scheme
%p: power of convolution term
%n: size of the set of masses
%iter_num: number of iterations you would like to perform
%index: form of network.
%0- two points, trivial network with two stations
%1- network crossing at the center
%2- network of 7 points, hexagon trivlet pattern
%3- a consecutive line of stations
function Test_all(index, n, iter_num)
    switch index
        case 0
            Test_trivial_network(iter_num);
        case 1 
            Test_trivial_network1(n,iter_num);
        case 2
            Test_trivial_network2(n,iter_num);
        case 3
            Test_trivial_network3(n, iter_num);
        case 4
            Test_no_network(n, iter_num);
    end
end
function res = Test_no_network(n, iter_num)
    n = 2;
    p = 3;
    h = 0.005;
    alpha = 0.1;
    theta = 0.05;
    sigma = 0.05;
    X = [0.49 0.51; 0.51 0.49];
    M = ones(1,n)/n;
    Y = [0.48 0.5; 0.52 0.5];
    Adj = zeros(2);
    color_mat = rand(n,3);
    res = euler_wrap(X, Y, Adj, M, iter_num, p, h, alpha, theta, sigma, 0.02, color_mat);
end

function res = Test_trivial_network(iter_num)
    n = 2;
    p = 3;
    h = 0.005;
    alpha = 0.1;
    theta = 0.05;
    sigma = 0.05;
    X = [0.35 0.55; 0.65 0.55];
    M = ones(1,n)/n;
    Y = [0.3 0.5; 0.7 0.5];
    Adj = zeros(2);
    Adj(1,2) = 1;
    Adj = or(Adj, Adj');
    color_mat = rand(n,3);
    res = euler_wrap(X, Y, Adj, M, iter_num, p, h, alpha, theta, sigma, 0.02, color_mat);
end

function res = Test_trivial_network1(n,iter_num)
    p = 3;
    h = 0.005;
    alpha = 0.1;
    theta = 0.3;
    sigma = 0.05;
    M = ones(1,n)/n;
    X = rand(n,2);
    Y = [0.3 0.7;0.7 0.7;0.5 0.5;0.3 0.3; 0.7 0.3];
    Adj = zeros(5);
    Adj(1,3) = 1;
    Adj(2,3) = 1;
    Adj(3,4) = 1;
    Adj(3,5) = 1;
    Adj = or(Adj, Adj');
    color_mat = rand(n,3);
    res = euler_wrap(X, Y, Adj, M, iter_num, p, h, alpha, theta, sigma, 0.02,color_mat);
end


function res = Test_trivial_network2(n,iter_num)
    p = 3;
    h = 0.005;
    alpha = 0.1;
    theta = 0.03;
    sigma = 0.05;
    M = ones(1,n)/n;
    X = rand(n,2);
    Y = [0.5 0.7;0.3268 0.6;0.3268 0.4;0.5 0.3; 0.6732 0.4;0.6732 0.6; 0.5 0.5];
    Adj = zeros(size(Y));
    Adj(1,2) = 1;
    Adj(2,3) = 1;
    Adj(3,4) = 1;
    Adj(4,5) = 1;
    Adj(5,6) = 1;
    Adj(6,1) = 1;
    Adj(1,7) = 1;
    Adj(3,7) = 1;
    Adj(5,7) = 1;
    Adj = or(Adj, Adj');
    color_mat = rand(n,3);
    res = euler_wrap(X, Y, Adj, M, iter_num, p, h, alpha, theta, sigma, 0.02, color_mat);
end

function res = Test_trivial_network3(n,iter_num)
    p = 3;
    h = 0.005;
    alpha = 0.2;
    theta = 0.002;
    sigma = 0.05;
    M = ones(1,n)/n;
    X = rand(n,2);
    Y = 5:13;
    Y = [Y', ones(9,1)*10]/20; 
    Adj = zeros(size(Y));
    Adj(1,2) = 1;
    Adj(2,3) = 1;
    Adj(3,4) = 1;
    Adj(4,5) = 1;
    Adj(5,6) = 1;
    Adj(6,7) = 1;
    Adj(7,8) = 1;
    Adj(8,9) = 1;
    Adj = or(Adj, Adj');
    color_mat = rand(n,3);
    res = euler_wrap(X, Y, Adj, M, iter_num, p, h, alpha, theta, sigma, 0.02, color_mat);
end