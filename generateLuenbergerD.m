% generateLuenbergerD   Generates the state matrix D to use in a Luenberger
% observer given the desired eigenvalues.
%
%   D = generateLuenbergerD(eigenvals) returns a real m X m matrix D
%   representing the state matrix of a Luenberger observer. 'eigenvals' is
%   a complex m X 1 matrix containing the desired eigenvalues of the
%   observer. It is assumed all complex eigenvalues are in conjugate pairs.

function D = generateLuenbergerD(eigenvals)
% Written by Nils Wilhelmsen, October 2020
%
% Function description: Given an array of eigenvalues 'eig', function
% returns real matrix D to use in Luenberger observer.
%
% Function presumption: It is assumed all eigenvalues in 'eig' are unique, 
% and that complex eigenvalues always come in conjugate pairs.

%% Step 1: Set up seperate arrays for complex and real eigenvalues
[eig_complex, eig_real] = splitEigsComplexReal(eigenvals);

%% Step 2: Sort complex eigenvalue array into conjugate pairs
if(~isnan(eig_complex))
    eig_complex = sort(eig_complex);
end
%% Step 3: Set up cell array containing 2X2 matrices for each complex conjugate pair, and scalars for each real eigenvalue
eig_cell = generateCellEigComplexReal(eig_complex, eig_real);

%% Step 4: Set up D as real bulk diagonal matrix 
D = blkdiag(eig_cell{:});
end
