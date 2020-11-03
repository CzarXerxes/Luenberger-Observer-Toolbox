% computeLFromLuenberger    Computes the observer gain corresponding to a
% given observer input matrix F and left-inverse transformation T_star.
%
%   L = computeLFromLuenberger(T_star, F) returns n X 1 real matrix L
%   representing the observer gain for an output-error injection observer, 
%   computed from the n X m real matrix T_star, representing the linear 
%   left-inverse Luenberger observer, and m X 1 real matrix F, representing 
%   the input matrix of the Luenberger observer.

function L = computeLFromLuenberger(T_star, F)
% Written by Nils Wilhelmsen, October 2020
%
% Function description:
%   - This function computes the observer gain L corresponding to
%   Luenberger observer with input matrix F and left-inverse Luenberger
%   observer T_star.
%
% Function presumption:
%   - The input arguments are both real matrices of appropriate dimension.

%% Step 1: Compute L
L = T_star*F;
end