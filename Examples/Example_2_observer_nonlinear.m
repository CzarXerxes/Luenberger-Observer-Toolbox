%% Clean up workspace before starting
clear all
close all
%% System definition
% System model
f = @(x) [x(2).^3; -x(1);];
g = @(x) [x(1) - x(1); x(2) - x(2);];

% System measurement
h = @(x) x(1);

% System test input
u0 = @(t) 0*t;

% System model dimension
dimx = 2;

%% Luenberger observer design
[b,a] = besself(3, 2*pi);
eigen = roots(a);

D = generateLuenbergerD([eigen(1) eigen(2) eigen(3)]);
F = [1;1;1;];

%% Simulation parameters
% Number of simulations
nsims = 10;
% Length of simulations
tsim = 40;
% Length of time step to store results as
dt = 1e-2;
% Array of initial conditions; 
w0_array = zeros(dimx + size(D,2),nsims);
for idx=1:nsims
    w0_array(1,idx) = 0.1*idx;
end

%% Estimate T_star

% Specify neural network parameters
hidden_layers_size_arr = [25 25];
bias_connect = [1; 0; 0;];

% Obtain T_star
T_star_net = computeNonlinearLuenbergerTnn(f, g, h, dimx, u0, D, F, w0_array, nsims, tsim, dt, hidden_layers_size_arr, bias_connect);


%% Test observer
w0_test = [rand(1) - 1; rand(1) - 1; 0; 0; 0;];
[tq_test, w_test] = performMultipleLuenbergerSimulations(f,g,h,dimx,u0,D,F,w0_test,1,tsim,dt);

x_hat = zeros(length(tq_test),2);
for jdx=1:length(tq_test)
    x_hat(jdx,:) = T_star_net(w_test(jdx,3:5)');
end

%% Plots
figure
hold on
plot(tq_test, x_hat(:,1));
plot(tq_test, w_test(:,1));

figure
hold on
plot(tq_test, x_hat(:,2));
plot(tq_test, w_test(:,2));
