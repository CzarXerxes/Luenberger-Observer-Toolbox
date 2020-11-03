% estimateLinearXOutputError    Estimates the state x of a SISO LTI ODE
% system via an output-error injection observer.
%
%   [t, x_hat] = estimateLinearXOutputError(A,B,C,L,x0,tsim,u) calculates
%   the n-dimensional state estimate x_hat at times t of the plant with
%   real n X n state matrix A, real n X 1 input matrix B, real 1 X n
%   output matrix C and scalar input u, through simulating an output-error
%   injection observer with gain L together with plant. The plant is
%   initialized at x0, wheres the observer is initialized at the origin.
%   The simulation is performed from t=0 to t=sim.
%
%   [t, x_hat, x] = estimateLinearXOutputError(A,B,C,L,x0,tsim,u)
%   calculates the n-dimensional state estimate x_hat at times t, and the
%   corresponding plant state x. The input arguments are as specified
%   above. 

function [t, x_hat, x] = estimateLinearXOutputError(A,B,C,L,x0,tsim,u)
% Written by Nils Wilhelmsen, October 2020
%
% Function description:
%   - Given the SISO LTI ODE system with state X, input U, output Y:
%
%               X'(t) = AX(t) + BU(t)
%                Y(t) = CX(t) + HU(t)
%
%   - Output-error injection observer implemented by copying state and
%   adding term consisting of observer gain L multiplied by
%   output-estimation error:
%
%             X_hat'(t) = AX_hat(t) + BU(t) + L[Y(t) - CX_hat(t) - HU(t)]
%
%   - Function simulates plant and output-error injection observer
%   simulatenously, giving the option of returning both X and X_hat.
%
% Function presumption:
%
%   - It is assumed the first four input arguments are matrices of
%   appropriate dimension, i.e. A should be n X n, B should be n X 1, C
%   sould be 1 X n, and L should be n X 1. Here n is the dimension of X.
%   For the observer to produce convergent estimates, it must be ensured
%   that A-LC is Hurwitz.
%
%   - The fifth input argument should be an n X 1 matrix representing the
%   plant initial condition.
%
%   - The sixth input argument tsim should be a positive real number
%   representing the duration of the simulation.
%
%   -The seventh input argument u should be a function handle representing
%   the plant input signal.
%
%   - The output argument t is a q X 1 matrix containing the sampling
%   points during the simulation, where q is the number of sampling points
%   taken. x_hat is a q X n matrix containing the state estimates for each
%   sampling point, and x is a q X n matrix containing the plant state at
%   each sampling point.

%% Step 1: Set up combined plant and observer matrices
A_w = [A zeros(size(A,1)); L*C (A-L*C)];
B_w = [B; B;];

%% Step 2: Set up simulation data structures
tspan = [0 tsim];
w0 = [x0; zeros(size(A,1),1);];

%% Step 3: Perform integration
[t,w] = ode45(@(t,w) A_w*w + B_w*u(t), tspan, w0);

%% Step 4: Distribute results into respective storage arrays
x = w(1:end, 1:size(A,1));
x_hat = w(1:end, (size(A,1)+1):(2*size(A,1)));

end