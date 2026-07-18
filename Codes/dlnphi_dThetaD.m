function dlnphi_dThetaD = dlnphi_dThetaD(Z, A, B, Bi, Aij_matrix, ...
                                                         dZ_dThetaD, dA_dThetaD, dB_dThetaD, dAij_dThetaD_all, ...
                                                         u1, u2)
% COMPUTE_DLNPHI_DTHETAD_MATRIX computes ∂ln(φ_i)/∂Θ_Dk for all i and k
%
% Inputs:
%   Z                   - scalar compressibility factor
%   A, B                - scalar EOS parameters
%   Bi                  - Nx1 vector of component-specific B_i
%   Aij_matrix          - Nx1 vector of A_ij (diagonal terms A_ii or mean values)
%   dZ_dThetaD          - 6x1 vector of ∂Z/∂Θ_Dk
%   dA_dThetaD          - 6x1 vector of ∂A/∂Θ_Dk
%   dB_dThetaD          - 6x1 vector of ∂B/∂Θ_Dk
%   dAij_dThetaD_all    - Nx6 matrix: each row is ∂A_ij/∂Θ_Dk for component i
%   u1, u2              - EOS constants (e.g., PR: 1±sqrt(2))
%
% Output:
%   dlnphi_dThetaD      - Nx6 matrix of ∂ln(φ_i)/∂Θ_Dk

    N = length(Bi);
    dlnphi_dThetaD = zeros(N, 6);

    Zpu1B = Z + u1 * B;
    Zpu2B = Z + u2 * B;
    ln_ratio = log(Zpu1B / Zpu2B);

    for i = 1:N
        for k = 1:6
            dZ = dZ_dThetaD(k);
            dA = dA_dThetaD(k);
            dB = dB_dThetaD(k);
            dAij = dAij_dThetaD_all(i, k);

            Bi_i = Bi(i);
            Aij = Aij_matrix(i);

            % --- Term 1: B_i/B and Z ---
            term1 = (Bi_i / B) * dZ - (Bi_i * (Z - 1) / B^2) * dB;

            % --- Term 2: ln(Z - B) derivative ---
            term2 = - (1 / (Z - B)) * (dZ - dB);

            % --- Term 3: nonlinear EOS term ---
            bracket1 = 2 * (Aij / A) - (Bi_i / B);
            bracket2 = (1 / Zpu1B) * (dZ + u1 * dB) - (1 / Zpu2B) * (dZ + u2 * dB);
            term3 = - (A / ((u1 - u2) * B)) * bracket1 * bracket2;

            % --- Term 4: derivative of the ln(...) term ---
            bracket31 = (2 / ((u1 - u2) * B)) *dAij;
            bracket32 =-2*Aij/((u1-u2)*A*B)*dA;
            bracket33 =A*Bi_i/((u1-u2)*B^3)*dB;
            term4 = - ln_ratio * (bracket31+bracket32+bracket33);
            % --- Term 5: derivative of the ln(...) term ---
            bracket51 = 2*Aij/A-Bi_i/B;
            bracket52 =(1/((u1-u2)*B)*dA - A/((u1-u2)*B^2)*dB);
            term5 = - ln_ratio * (bracket51*bracket52);

            % --- Total derivative ---
            dlnphi_dThetaD(i, k) = term1 + term2 + term3 + term4+term5;
        end
    end
end
