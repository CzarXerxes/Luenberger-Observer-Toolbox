% splitEigsComplexReal  Splits an array of eigenvalues containing possibly
% both complex and real values into one array containing only complex
% eigenvalues, and one array containing only real eigenvalues.
%
%   [eig_complex, eig_real] = splitEigsComplexReal(eigenvals) returns a
%   complex m_c X 1 matrix 'eig_complex' containing all the complex 
%   eigenvalues in 'eigenvals', and a real m_r X 1 matrix 'eig_real'
%   containing all the real eigenvalues. 'eigenvals' is an m X 1 complex
%   matrix, where it is assumed all complex eigenvalues come in conjugate
%   pairs. It is assumed m_c is an even number and m = m_c + m_r.

function [eig_complex, eig_real] = splitEigsComplexReal(eigenvals)
% Written by Nils Wilhelmsen, October 2020
%
% Function description: Given array of eigenvalues 'eigenvals', the function
% splits it into two arrays, 'eig_complex', containing only complex values,
% and 'eig_real', containing only real values. If 'eigenvals' only contains
% real (resp. complex) values, then 'eig_complex' (resp. 'eig_real')
% returns as NaN.
%
% Function presumption: The data structure 'eigenvals' is an array (1D
% matrix) containing real and/or complex eigenvalues. It contains at least
% one value.
%% Step 1: Declare complex and real eigenvalue arrays

%% Step 1: Declare eigenvalue data structures.
eig_complex = nan;
eig_real = nan;

%% Step 2: Initialize iteration indices
ind_complex = 1;
ind_real = 1;

%% Step 3: Iterate through 'eigenvals' array and split into 'eig_complex' and 'eig_real'
for val=eigenvals
    if imag(val) == 0
        eig_real(1,ind_real) = val;
        ind_real = ind_real + 1;
    else
        eig_complex(1,ind_complex) = val;
        ind_complex = ind_complex + 1;
    end
end
end