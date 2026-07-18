function dA_dThetaD = dA_dThetaD_vector(Theta_D)
% DA_DTHETAD_VECTOR computes ∂A/∂Θ_Dk as a 6x1 vector
%
% Input:
%   Theta_D - 6x1 vector [Θ_D1; Θ_D2; ...; Θ_D6]
%
% Output:
%   dA_dThetaD - 6x1 vector of partial derivatives

    dA_dThetaD = zeros(1, 6);
    dA_dThetaD(1) = 2 * Theta_D(1);
    dA_dThetaD(2) = 4 * Theta_D(2);
    dA_dThetaD(3) = -2 * Theta_D(4);
    dA_dThetaD(4) = -2 * Theta_D(3);
    % dA_dThetaD(5) = 0;
    % dA_dThetaD(6) = 0;
end
