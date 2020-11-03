% generateCellEigComplexReal    Generates a cell array containing 2X2 real
% matrices for each pair of complex conjugate eigenvalues passed in the
% arguments, and real scalar for each real eigenvalue.
%
%   eig_cell = generateCellEigComplexReal(eig_complex, eig_real) returns an
%   (m_c/2 + m_r)-dimensional cell array containing a real 2X2 matrix for
%   each complex conjugate pair in the complex m_c X 1 matrix 'eig_complex'
%   containing only complex conjugate pairs, and a real scalar for each
%   eigenvalue in the real m_r X 1 matrix 'eig_real'.

function eig_cell = generateCellEigComplexReal(eig_complex, eig_real)
% Written by Nils Wilhelmsen, October 2020
%
% Function description: Takes in a sorted array of complex eigenvalues,
% 'eig_complex', and an array of real eigenvalues, 'eig_real', and outputs
% a cell array containing 2X2 matrices for the complex eigenvalues and
% scalars for the real eigenvalues.
% 
% Function presumption: 'eig_complex' should only contain conjugate complex
% pairs and be sorted by these pairs. 'eig_real' should only contain real
% values. At least one of these arrays are not equal to NaN.

%% Step 1: Compute number of complex conjugate pairs and real eigenvalues
if(~isnan(eig_complex))
    num_complex_pairs = length(eig_complex)/2;
else
    num_complex_pairs = 0;
end

if(~isnan(eig_real))
    num_reals = length(eig_real);
else
    num_reals = 0;
end

%% Step 2: Declare cell array
eig_cell = cell(num_complex_pairs + num_reals,1);

%% Step 3: Set up 2X2 matrices for complex conjugate pairs and store in cell array
if(~isnan(eig_complex))
    temp = zeros(2);
    for h=1:2:2*num_complex_pairs
        % Set up matrix
        temp(1,1) = real(eig_complex(1,h));
        temp(1,2) = imag(eig_complex(1,h));
        temp(2,1) = imag(eig_complex(1,h+1));
        temp(2,2) = real(eig_complex(1,h+1));

        % Store matrix in cell array
        j =  (h + 1) / 2;
        eig_cell{j} = temp;
    end
end

%% Step 4: Store real eigenvalues in cell array
if(~isnan(eig_real))
    for j=(num_complex_pairs+1):(num_complex_pairs+num_reals)
        eig_cell{j} = eig_real(j-num_complex_pairs);
    end
end
