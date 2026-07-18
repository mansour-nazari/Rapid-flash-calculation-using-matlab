function Jac = compute_main_jacobian(Ai, h, g, Bi, dx_dTheta, dy_dTheta)
% COMPUTE_MAIN_JACOBIAN returns the 6x6 Jacobian of residuals ∂R_Dl/∂Θ_Dk
%
% Inputs:
%   Ai         - N×1 vector of A_i
%   h, g       - N×1 vectors of h_i and g_i
%   Bi         - N×1 vector of B_i
%   dx_dTheta  - N×6 matrix of ∂x_i/∂Θ_Dk^L
%   dy_dTheta  - N×6 matrix of ∂y_i/∂Θ_Dk^L
%
% Output:
%   Jac        - 6×6 Jacobian matrix

    N = length(Ai);
    Jac = zeros(6, 6);

    sqrt_A = sqrt(Ai);
    w1 = sqrt_A;
    w2 = sqrt_A .* h .* g;
    w3 = sqrt_A .* h.^2 .* g;
    w4 = sqrt_A .* g;
    w5 = Bi;

    W = [w1, w2, w3, w4, w5];  % N×5 weights

    % Rows 1 to 5
    for l = 1:5
        for k = 1:6
            Jac(l, k) = double(l == k) - sum(W(:, l) .* dx_dTheta(:, k));
        end
    end

    % Row 6 (R_D6)
    for k = 1:6
        Jac(6, k) = sum(dx_dTheta(:, k) - dy_dTheta(:, k));
    end
end
