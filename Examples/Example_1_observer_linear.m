%% Clean up the workspace before starting
clear all
close all
%% System definition

% System model
A = [-1 -2; 0 -3;];
B = [1; 2;];
C = [1 0];
H = 0;

% Specify initial conditions
x0 = [1; 3;];

% Input function
u = @(t) sin(t);

% Length of simulation
tsim = 10;

%% Luenberger observer
% Design observer
D = generateLuenbergerD([(-1 - 1i) (-1 + 1i) -2 ]);
F = [1; 1; 1;];

% Compute Luenberger forward and left-inverse transformations
[T, T_star] = computeLinearLuenbergerT(A,C,D,F);

% Simulate system with Luenberger observer
[t_luen,x_hat_luen,x_luen,z] = estimateLinearXLuenberger(A,B,C,H,D,F,x0,tsim,u,T,T_star);

% Luenberger observer plots
figure
hold on
plot(t_luen,x_luen(:,1));
plot(t_luen,x_hat_luen(:,1));

figure
hold on
plot(t_luen,x_luen(:,2));
plot(t_luen,x_hat_luen(:,2));

figure
hold on
plot(t_luen,z(:,1));
plot(t_luen,z(:,2));
plot(t_luen,z(:,3));
title('Observer states')

%% Output-error injection observer
% Obtain corresponding L matrix
L = computeLFromLuenberger(T_star, F);

% Simulate system with output-error injection observer
[t_oe, x_hat_oe, x_oe] = estimateLinearXOutputError(A,B,C,L,x0,tsim,u); 

% Output-error injection plots
figure
hold on
plot(t_oe,x_oe(:,1));
plot(t_oe,x_hat_oe(:,1));

figure
hold on
plot(t_oe,x_oe(:,2));
plot(t_oe,x_hat_oe(:,2));

