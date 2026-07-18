function dZ_dThetaD = compute_dZ_dThetaD(Z, A, B, u1, u2, dA_dThetaD, dB_dThetaD)
% COMPUTE_DZ_DTHETAD calculates ∂Z/∂Θ_Dk as a 6x1 vector
%
% Inputs:
%   Z            - compressibility factor (scalar)
%   A, B         - EOS parameters (scalars)
%   u1, u2       - EOS constants (e.g., PR: u1 = 1 + sqrt(2), u2 = 1 - sqrt(2))
%   dA_dThetaD   - ∂A/∂Θ_Dk, vector (6x1)
%   dB_dThetaD   - ∂B/∂Θ_Dk, vector (6x1)
%
% Output:
%   dZ_dThetaD   - ∂Z/∂Θ_Dk, vector (6x1)

    % Precompute intermediate terms
    a2 = (u1 + u2 - 1) * B-1;
    a1 = A - (u1 + u2) * B -(u1+u2- u1 * u2 )* B^2;
    a0 = -A * B - u1 * u2 * B^2 * (1 + B);

    % Precompute common prefactor
    denom = 3 * Z^2 + 2 * a2 * Z + a1;

    % Allocate output
    dZ_dThetaD = zeros(1, 6);

    for k = 1:6
        % Derivatives of a2, a1, a0
        da2 = (u1 + u2 - 1) * dB_dThetaD(k);

        da1 = dA_dThetaD(k) ...
            - u2 * dB_dThetaD(k) ...
            - 2 * (u1 + u2 - u1 * u2) * B * dB_dThetaD(k);

        da0 = -A * dB_dThetaD(k) ...
            - B * dA_dThetaD(k) ...
            - 2 * u1 * u2 * (1 + B) * B * dB_dThetaD(k) ...
            - u1 * u2 * B^2 * dB_dThetaD(k);

        % Chain rule for dZ/dThetaD_k
        numerator = Z^2 * da2 + Z * da1 + da0;
        dZ_dThetaD(k) = - numerator / denom;
    end
end
