function Jac = full_residual_jacobian(ThetaL, ThetaV, Z, A, B, a, alpha, b, h, g, z, K, Rgas, T, P, u1, u2, Aij)
% FULL_RESIDUAL_JACOBIAN computes the 6x6 Jacobian matrix ∂R_Dl/∂Θ_Dk^L
% using all internal helper functions
%
% INPUTS:
%   ThetaL     - 6x1 vector of Θ^L_Dk (liquid phase dimensionless Thetas)
%   ThetaV     - 6x1 vector of Θ^V_Dk (vapor phase dimensionless Thetas)
%   Z          - scalar compressibility factor
%   A, B       - EOS parameters (scalars)
%   a, alpha   - EOS component vectors (Nx1)
%   b, h, g    - component property vectors (Nx1)
%   z          - feed composition vector (Nx1)
%   K          - equilibrium ratios (Nx1)
%   Rgas, T, P - gas constant, temperature, pressure
%   u1, u2     - EOS constants (e.g., PR: 1±sqrt(2))
%   Aij        - Nx1 vector of A_ij values for each component
%
% OUTPUT:
%   Jac        - 6×6 Jacobian matrix of residuals

    N = length(a);

    % STEP 1: EOS parameter derivatives
    [Ai, Bi, dAij_dThetaD_all] = Jac_parametrs(a, alpha, b, h, g, P, T, Rgas);

    % STEP 2: Derivatives of A and B
    Theta_D = compute_dimensionless_Theta(ThetaL, P, T, Rgas);
    dA_dThetaD = dA_dThetaD_vector(Theta_D);
    dB_dThetaD = zeros(6,1); dB_dThetaD(5) = 1;

    % STEP 3: dZ/dThetaD
    dZ_dThetaD = compute_dZ_dThetaD(Z, A, B, u1, u2, dA_dThetaD, dB_dThetaD);

    % STEP 4: ln(phi) derivatives
    dlnphi = dlnphi_dThetaD(Z, A, B, Bi, Aij, dZ_dThetaD, dA_dThetaD, dB_dThetaD, dAij_dThetaD_all, u1, u2);

    % STEP 5: dThetaV/dThetaL
    dThetaV_dThetaL = dThetaDv_dThetaDl_matrix(ThetaL, ThetaV);

    % STEP 6: ∂lnφ^V/∂Θ^L using chain rule
    dlnphiV_dThetaL = zeros(N, 6);
    for i = 1:N
        dlnphiV_dThetaL(i,:) = dlnphi(i,:) * dThetaV_dThetaL;
    end

    % STEP 7: dK/dTheta
    dK_dTheta = compute_dK_dThetaDl_matrix(K, dlnphi, dThetaV_dThetaL, dlnphiV_dThetaL);

    % STEP 8: dx_i/dTheta and dy_i/dTheta
    theta6 = ThetaL(6);
    [dx_dTheta, dy_dTheta] = compute_dx_dy_dTheta(K, z, dK_dTheta, theta6);

    % STEP 9: Residual Jacobian
    Jac = compute_main_jacobian(Ai, h, g, Bi, dx_dTheta, dy_dTheta);
end
