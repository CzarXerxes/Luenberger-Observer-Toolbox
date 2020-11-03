% computeLinearLuenbergerG  Generates the G matrix necessary to implement
% Luenberger observers for SISO LTI ODE systems with an input signal.
%
%   G = computeLinearLuenbergerG(B,H,F,T) returns an m X 1 real matrix G
%   representing the matrix to be multiplied by the input signal U in a
%   linear Luenberger observer. Here B is an n X 1 real matrix
%   representing the input matrix of the plant, H is a real scalar
%   representing the transfer constant of the plant, F is an m X 1 real
%   matrix representing the input matrix of the observer which multiplies
%   the output signal Y, and T is an m X n real matrix representing the
%   forward linear Luenberger transformation for the given plant and
%   Luenberger observer.

function G = computeLinearLuenbergerG(B,H,F,T)
% Written by Nils Wilhelmsen, October 2020
%
% Function description: Computes the G necessary to implement Luenberger
% observer for linear ODEs with input signal.
% 
% Function presumption: All matrices B,H,F,T are real matrices of suitable
% size
%
%% Step 1: Compute G
G = T*B - F*H;
end

