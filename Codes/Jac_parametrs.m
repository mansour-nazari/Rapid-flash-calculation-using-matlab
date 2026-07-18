function [Ai, Bi, dAij_dTheta] = Jac_parametrs(a, alpha, b, h, g, P, T, R)
% JAC_PARAMETRS computes EOS parameters Ai, Bi and their sensitivity to Θᴸ
%
% Inputs:
%   a      - vector of a_i
%   alpha  - vector of alpha_i
%   b      - vector of b_i
%   h, g   - vectors of h_i and g_i
%   P      - pressure (scalar)
%   T      - temperature (scalar)
%   R      - gas constant (scalar)
%
% Outputs:
%   Ai            - vector of A_i values
%   Bi            - vector of B_i values
%   dAij_dTheta   - N×6 matrix of ∂A_ij/∂Θ_kᴸ for i = j

    % Compute A_i and B_i
    Ai = (P / (R * T)^2) * (a .* alpha);
    Bi = (P / (R * T))   * b;

    % Sensitivity matrix ∂A_ij/∂Θ_kᴸ for i = j
    N = length(a);
    dAij_dTheta = zeros(N, 6);
    
    for i = 1:N
        Aij = Ai(i);  % A_ii
        sqrt_Aij = sqrt(Aij);
        dAij_dTheta(i, :) = [ ...
            sqrt_Aij, ...
            2 * sqrt_Aij * h(i) * g(i), ...
            -sqrt_Aij * g(i), ...
            -sqrt_Aij * h(i)^2 * g(i), ...
            0, ...
            0 ...
        ];
    end
end
