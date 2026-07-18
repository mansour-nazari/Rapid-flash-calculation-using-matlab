function [dx_dTheta, dy_dTheta] = compute_dx_dy_dTheta(K, z, dK_dTheta, theta6)
% COMPUTE_DX_DY_DTHETA calculates ∂x_i/∂Θ_Dk^L and ∂y_i/∂Θ_Dk^L for all i, k
%
% Inputs:
%   K          - N×1 vector of K_i
%   z          - N×1 vector of z_i
%   dK_dTheta  - N×6 matrix of ∂K_i/∂Θ_Dk^L
%   theta6     - scalar, Θ^L_D6 (liquid molar fraction)
%
% Outputs:
%   dx_dTheta  - N×6 matrix of ∂x_i/∂Θ_Dk^L
%   dy_dTheta  - N×6 matrix of ∂y_i/∂Θ_Dk^L

    N = length(K);
    dx_dTheta = zeros(N, 6);
    dy_dTheta = zeros(N, 6);

    for i = 1:N
        for k = 1:6
            delta6k = double(k == 6);
            denom = (theta6 + K(i) * (1 - theta6))^2;

            % dx_i/dTheta_k^L
            num_x = (1 - K(i)) * delta6k + (1 - theta6) * dK_dTheta(i, k);
            dx_dTheta(i, k) = -z(i) * num_x / denom;

            % dy_i/dTheta_k^L
            num_y = K(i) * (1 - K(i)) * delta6k + theta6 * dK_dTheta(i, k);
            dy_dTheta(i, k) = z(i) * num_y / denom;
        end
    end
end
