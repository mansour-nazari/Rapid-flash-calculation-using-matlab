function [Psi, a_m, b_m] = compute_psi_am_bm_from_theta(a, alpha, h, g, Theta)
% COMPUTE_PSI_AM_BM_FROM_THETA computes Ψ_i, a_m, b_m from Θⁱᴸ values
%
% Inputs:
%   a, alpha, h, g - column vectors (length N)
%   Theta - vector of length 6 containing Θ₁ᴸ to Θ₆ᴸ
%
% Outputs:
%   Psi - vector of Ψᵢ values (N x 1)
%   a_m - scalar
%   b_m - scalar

    N = length(a);
    Psi = zeros(N, 1);
    
    for i = 1:N
        sqrt_term = sqrt(a(i) * alpha(i));
        Psi(i) = sqrt_term * Theta(1) + ...
                 2 * sqrt_term * h(i) * g(i) * Theta(2) - ...
                 sqrt_term * g(i) * Theta(3) - ...
                 sqrt_term * h(i)^2 * g(i) * Theta(4);
    end

    a_m = Theta(1)^2 + 2 * Theta(2)^2 - 2 * Theta(3) * Theta(4);
    b_m = Theta(5);
end
