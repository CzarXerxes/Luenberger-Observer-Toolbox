% computeLinearLuenbergerT  Compute the linear Luenberger transformation
% of a SISO LTI ODE system, and the corresponding left-inverse 
% transformation, whenever these exist.
%
%   T = computeLinearLuenbergerT(A,C,D,F) returns a real m X n matrix T 
%   representing the forward linear Luenberger transformation mapping the 
%   n-dimensional plant state X, with corresponding state matrix A and 
%   output matrix C, into the m-dimensional target system state Z, with
%   corresponding state matrix D and input matrix F, whenver it exists. 
%   If A and D have any common eigenvalues then T does not exist and NaN is
%   returned instead.
%
%   [T, T_star] = computeLinearLuenbergerT(A,C,D,F) returns a real m X n
%   matrix T representing the forward linear Luenberger transformation as
%   described above, and a real n X m matrix T_star representing the
%   corresponding left-inverse transformation, whenever these exist. If
%   either (A,C) is not observable or (D,F) is not controllable, T_star 
%   does not exist and NaN is returned instead.

function [T,T_star] = computeLinearLuenbergerT(A,C,D,F)
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
%   - The function computes the linear mapping T connecting the state Z and 
%      X via:
%
%               Z(t) = TX(t)
%
%     and the corresponding left-inverse T_star, if these exist.
%
% Function presumption:
%
%   - It is assumed the input arguments are real matrices of appropriate
%     dimension, i.e. A should be n X n, C should be 1 X n, D should be m X
%     m, and F should be m X 1. Here n is the dimension of X, and m is the
%     dimension of Z.
%
%   - The forward transformation T only exists if A and D have no common
%     eigenvalues. If the matrices A and D have common eigenvalue(s), NaN 
%     is returned for both T and T_star.
%
%   - Given that T exists, the left-inverse T_star only exists if the pair
%     (A,C) is observable and (D,F) is controllable. If either of these two 
%     criteria are not satisfied, NaN is returned for T_star.

%% Step 1: Check whether T exists (A,D have no common eigenvalues)
eig_A = eig(A);
eig_D = eig(D);
common_eigs = ~isempty(intersect(eig_A,eig_D));
if(common_eigs)
    T = nan;
    T_star = nan;
    return;
end
%% Step 2: T exists, compute it
T = sylvester(-D,A,F*C);

%% Step 3: Check whether T_star exists( (A,C) observable and (D,F) controllable)
Ob = obsv(A,C);
unob = length(A) - rank(Ob);

Co = ctrb(D,F);
unco = length(D) - rank(Co);

if(unob + unco ~= 0)
    T_star = nan;
    return;
end
%% Step 4: Compute inverse
T_star = pinv(T);
end