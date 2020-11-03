% computeNonlinearLuenbergerTnn  Numerically estimate the
% nonlinear Luenberger transformation of a SISO input-affine nonlinear
% system with static transformation, and the corresponding left-inverse.
%
%   T_star_net =
%   computeNonlinearLuenbergerTnn(f,g,h,dimx,D,F,w0_array,nsims,tsim,dt,...
%   hidden_layers_size_arr,bias_connect) returns a neural network with 
%   m input nodes and n output nodes, representing the left-inverse 
%   nonlinear Luenberger transformation mapping the m-dimensional observer 
%   state Z, with corresponding state matrix D and input matrix F, into the 
%   n-dimensional plant state X, with state function f,input-affine input 
%   function g and output function h. 'u0' is a function handle
%   representing the test input to be used for training the neural network.
%   'dimx' denotes the dimension of X, i.e. n. 'w0_array' is an (n+m) X 
%   'nsims' matrix representing the initial conditions to generate 
%   simulation data for training the neural network(s). 'tsim' is a 
%   positive real number representing the time the simulations should go on 
%   for. 'dt' denotes the time-interval the simulation data should be 
%   sampled at after the simulation is performed, for training the neural 
%   network. 'hidden_layers_size_arr' is a 1 X p matrix containing positive 
%   integers representing the size of the hidden layers, where p is the 
%   number of hidden layers. 'bias_connect' is a (p+1) X 1 boolean matrix 
%   representing which layers should have a bias connected to it. Conjugate 
%   gradient backpropagation with Polak-Ribiére updates is used as the 
%   training function by default.
%
%   [T_star_net, T_net] =
%   computeNonlinearLuenbergerTnn(f,g,h,dimx,D,F,w0_array,nsims,tsim,dt,...
%   hidden_layers_size_arr,bias_connect) returns a neural network
%   'T_star_net' representing the left-inverse Luenberger transformation 
%   described above, and also a neural network 'T_net' with n input nodes 
%   and m output nodes representing the corresponding forward 
%   transformation. The input arguments are as described above

function [T_star_net, T_net] = computeNonlinearLuenbergerTnn(f,g,h,dimx, u0, D,F, w0_array, nsims, tsim, dt, hidden_layers_size_arr, bias_connect)
% Written by Nils Wilhelmsen, October 2020
%
% Function description: 
%   -Given the SISO input-affine nonlinear system with state X, input U,
%   output Y:
%
%       X'(t) = f(X(t)) + g(X(t))U(t)
%        Y(t) = h(X(t))
%
%   -Given target system with state Z, driven by Y and U:
%
%       Z'(t) = DZ(t) + FY(t) + \phi(Z)(U(t) - U_0(t))
%
%   where U_0 is a test input.
%
%   -The function numerically approximates the nonlinear left-inverse 
%   mapping T_star corresponding to the given U_0, connecting the 
%   state-estimate X_hat to Z via
%
%           X_hat = T_star(Z)
%
%   It can also return the corresponding forward transformation T.

%% Step 1: Generate simulation data
[tq, output_data] = performMultipleLuenbergerSimulations(f,g,h,dimx,u0,D,F, w0_array,nsims,tsim,dt);

%% Step 2: Preprocess the data
% Remove initial data due to noninjectivity of transformation
k=3;
t_c = k/min(abs(real(eig(D))));
I = find(tq < t_c);
idx = max(I);
output_data = output_data(idx:end,:,:);

regression_data = zeros(nsims*size(output_data,1),size(output_data,2));

for jdx=1:nsims
    regression_data( ((jdx-1)*size(output_data,1)+1):(jdx*size(output_data,1)), :) = output_data(:,:,jdx);
end

%% Step 3: Set up and train the neural network

T_star_net = feedforwardnet(hidden_layers_size_arr);
T_star_net.biasConnect = bias_connect;
T_star_net.trainFcn = 'traincgp';
T_star_net = train(T_star_net, regression_data(:,(dimx+1):(dimx+size(D,2)))' , regression_data(:, 1:dimx)');

if(nargout == 2)
    T_net = feedforwardnet(hidden_layers_size_arr);
    T_star_net.biasConnect = bias_connect;
    T_star_net.trainFcn = 'traincgp';
    T_net = train(T_net,regression_data(:, 1:dimx)'  , regression_data(:,(dimx+1):(dimx+size(D,2)))' );
end
end
