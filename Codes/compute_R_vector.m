function R = compute_R_vector(a, alpha, h, g, b, x, Theta, z, K)
% COMPUTE_R_VECTOR computes the residual vector R (6x1)
%
% Inputs:
%   a, alpha, h, g, b, x - column vectors (Nx1)
%   Theta - vector of Theta_1^L to Theta_6^L (6x1)
%   z     - feed composition vector z_i (Nx1)
%   K     - equilibrium ratios K_i (Nx1)
%
% Output:
%   R - residual vector [R1; R2; ...; R6] (6x1)

    N = length(x);
    
    % Precompute terms
    sqrt_a_alpha = sqrt(a .* alpha);           % √(a_i * alpha_i)
    term1 = sum(x .* sqrt_a_alpha);
    term2 = sum(x .* sqrt_a_alpha .* h .* g);
    term3 = sum(x .* sqrt_a_alpha .* h.^2 .* g);
    term4 = sum(x .* sqrt_a_alpha .* g);
    term5 = sum(x .* b);

    % Residuals R1 to R5
    R = zeros(1,6);
    R(1) = Theta(1) - term1;
    R(2) = Theta(2) - term2;
    R(3) = Theta(3) - term3;
    R(4) = Theta(4) - term4;
    R(5) = Theta(5) - term5;

    % R6: summation over i of (z_i * (1 - K_i) / (K_i + Theta6 * (1 - K_i)))
    Theta6 = Theta(6);
    denom = K + Theta6 .* (1 - K);
    R(6) = sum(z .* (1 - K) ./ denom);

end
