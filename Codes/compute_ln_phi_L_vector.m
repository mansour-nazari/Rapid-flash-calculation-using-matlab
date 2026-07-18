function ln_phi = compute_ln_phi_L_vector(Z, A, B, Bi, Aij, u1, u2)
% COMPUTE_LN_PHI_L_VECTOR calculates ln(phi_i^L) for all components
%
% Inputs:
%   Z    - scalar compressibility factor in liquid phase
%   A    - scalar EOS parameter A
%   B    - scalar EOS parameter B
%   Bi   - vector of component-specific B_i (Nx1)
%   Aij  - matrix of interaction parameters A_ij (NxN)
%   u1, u2 - EOS constants (e.g., PR: 1±sqrt(2))
%
% Output:
%   ln_phi - vector of ln(phi_i^L) values (Nx6)

    N = length(Bi);
    ln_phi = zeros(N,1);

    % Precompute log term
    log_term = log((Z + u1 * B) / (Z + u2 * B));
    coeff = A / ((u1 - u2) * B);

    for i = 1:N
        sum_Aij = Aij(i,:);  % Sum over j for fixed i (∑_j A_ij)
        term1 = (Bi(i) / B) * (Z - 1);
        term2 = log(Z - B);
        bracket = (2 * sum_Aij / A) - (Bi(i) / B);
        term3 = coeff * bracket * log_term;

        ln_phi(i) = term1 - term2 - term3;
    end
end
