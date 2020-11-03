% performMultipleLuenbergerSimulations    Runs and outputs the results from 
% multiple simulations of an input-affine nonlinear system driving a 
% Luenberger observer target system.
%
%   [tq, output_data] =
%   performMultipleLuenbergerSimulations(f,g,h,dimx,u,D,F,w0_array,...
%   nsims,tsim,dt) returns the output of 'nsims' simulations lasting for
%   time 'tsim' of the n-dimensional state affine nonlinear system with 
%   state function 'f', input function 'g', output function 'g' and input 
%   'u', driving a Luenberger observer with m X m real state matrix D and 
%   m X 1 real input matrix F driven by output from plant. 'w0_array' is an
%   (n+m) X 'nsims' real matrix representing the initial conditions for the
%   plant and observer target system for the simulations to be performed.
%   'dt' is the time step the simulation results are resampled into before
%   getting returned.

function [tq, output_data] = performMultipleLuenbergerSimulations(f,g,h,dimx,u,D,F,w0_array,nsims,tsim,dt)
% Written by Nils Wilhelmsen, October 2020
%
% Function description: 
%


%% Step 1: Set up necessary data structures
tspan = [0 tsim];
tq = 0:dt:tsim;
output_data = zeros(length(tq), dimx + size(D,2), nsims);

%% Step 2: Perform simulations
for idx = 1:nsims
    w0 = w0_array(:,idx);
    [t,w] = ode45( @(t,w) nonlin(t,w,f,g,h,u,D,F), tspan, w0);
    wq = interp1(t,w,tq);
    output_data(:,:,idx) = wq;
end

%% Simulation function definition
    function dwdt = nonlin(t,w,f,g,h,u,D,F)
        x = w(1:dimx);
        dwdt = [f(x); 
                zeros(size(D,2),1);] + [zeros(dimx) zeros(dimx, size(D,2));
                                        zeros(size(D,2), dimx) D;]*w + [zeros(dimx,1);
                                                                                F;]*h(x) + [g(x);
                                                                                            zeros(size(D,2),1);]*u(t);                                                                     
    end
end