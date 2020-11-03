% estimateLinearXLuenberger     Estimates the state x of a SISO LTI ODE 
% system via a Luenberger observer. 
%
%   [t, x_hat] = estimateLinearXLuenberger(A,B,C,H,D,F,x0,tsim,u,T,T_star)
%   calculates the n-dimensional state estimate x_hat at times t of the
%   plant with real n X n state matrix A, real n X 1 input matrix B, 
%   real 1 X n output matrix C, scalar feedthrough constant H and scalar 
%   input signal u, through simulating an m-dimensional target system with 
%   state matrix D and input matrix F, driven by the plant output and input 
%   signals. The plant is initialized at x0, whereas the observer target 
%   system is initialized at the origin. The simulation is performed from
%   t=0 to t=tsim. T is a real m X n matrix representing the linear forward
%   Luenberger transformation corresponding to the specified plant and
%   target system. T_star is a real n X m matrix representing the
%   corresponding left-inverse transformation.
%
%   [t,x_hat,x] = estimateLinearXLuenberger(A,B,C,H,D,F,x0,tsim,u,T,T_star)
%   caclulates the n-dimensional state estimate x_hat at times t, and
%   corresponding plant state x. The input arguments are as specified
%   above.
%
%   [t,x_hat,x,z] =
%   estimateLinearXLuenberger(A,B,C,H,D,F,x0,tsim,u,T,T_star) calculates
%   the n-dimensional state estimate x_hat at times t, the corresponding
%   plant state x, and the m-dimensional observer state z. The input
%   arguments are as specified above.

function [t, x_hat, x, z] = estimateLinearXLuenberger(A,B,C,H,D,F,x0,tsim,u,T,T_star)
% Written by Nils Wilhelmsen, October 2020
%
% Function description:
%   - Given the SISO LTI ODE system with state X, input U, output Y:
%
%               X'(t) = AX(t) + BU(t)
%                Y(t) = CX(t) + HU(t)
%
%   - Given target sytem with state Z, driven by Y and U:
%
%               Z'(t) = DZ(t) + FY(t) + GU(t)
%
%   - The function computes the state esimate X_hat corresponding to the
%   plant state X via the left-inverse transformation:
%
%               X_hat = T_star * Z
%
% Function presumption:
%
%   - It is assumed the first six input arguments are real matrices of
%   appropriate dimension, i.e. A should be n X n, B should be n X 1, C
%   should be 1 X n, H should be 1 X 1, D should be m X m and F should be m
%   X 1. Here n is the dimension of X, and m is the dimension of Z.
%   
%   - The seventh input argument x0 should be an n X 1 matrix representing
%   the plant initial condition.
%
%   - The eighth input argument tsim should be a positive real number 
%   representing the duration of the simulation.
%
%   - The ninth input argument u should be a function handle representing the
%   plant input signal.
%
%   -The last two input argumnets should be real matrices of appropriate
%   dimension representing the Luenberger forward and left-inverse
%   transformation corresponding to the supplied plant and target system
%   matrices, i.e. T should be an m X n real matrix and T_star should be a
%   n X m real matrix.
%
%   - The output argument t is a q X 1 matrix containing the sampling 
%   points during the simulation, where q is the number of sampling points 
%   taken. x_hat is a q X n matrix containing the state estimates for each 
%   sampling point, x is a q X n matrix containing the plant state at each 
%   sampling point and z is a q X m matrix containing the observer state at 
%   each sampling point. 

%% Step 1: Compute G matrix
G = computeLinearLuenbergerG(B,H,F,T);

%% Step 2: Set up combined plant and observer matrices
A_w = [A zeros(size(A,1),size(D,2)); F*C D;];
B_w = [B; G+F*H;];

%% Step 3: Set up simulation data structures
tspan = [0 tsim];
w0 = [x0; zeros(size(D,2),1);];

%% Step 4: Perform integration
[t,w] = ode45(@(t,w) A_w*w + B_w*u(t), tspan, w0);

%% Step 5: Distribue results into respective storage arrays
x = w(1:end, 1:size(A,2));
z = w(1:end, (size(A,2)+1):(size(A,2) + size(D,2)));
x_hat = (T_star*z(:,1:end)')';

%% Future extension: Allow T, T_star to be automatically computed if not passed in as arguments
% [T,T_star] = computeLinearLuenbergerT(A,C,D,F);
% if(isnan(T(1,1)) || isnan(T_star(1,1)))
%     fprintf("Forward and/or backwards transformation does not exist!\n");
%     t=nan;
%     x=nan;
%     z=nan;
%     x_hat=nan;
%     return;
% end
end