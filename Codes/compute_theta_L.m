function Theta = compute_theta_L(a, alpha, h, g, x, b, L)
% COMPUTE_THETA_L computes an array Theta(1:6) based on the given formulas:
%
% Θ₁ᴸ = ∑ xᵢ sqrt(aᵢ αᵢ)
% Θ₂ᴸ = ∑ xᵢ sqrt(aᵢ αᵢ) hᵢ gᵢ
% Θ₃ᴸ = ∑ xᵢ sqrt(aᵢ αᵢ) hᵢ² gᵢ
% Θ₄ᴸ = ∑ xᵢ sqrt(aᵢ αᵢ) gᵢ
% Θ₅ᴸ = ∑ xᵢ bᵢ
% Θ₆ᴸ = L (input scalar)

    N = length(a);
    
    % Initialize Theta vector
    Theta = zeros(1, 6);
    
    for i = 1:N
        sqrt_a_alpha = sqrt(a(i) * alpha(i));
        Theta(1) = Theta(1) + x(i) * sqrt_a_alpha;
        Theta(2) = Theta(2) + x(i) * sqrt_a_alpha * h(i) * g(i);
        Theta(3) = Theta(3) + x(i) * sqrt_a_alpha * h(i)^2 * g(i);
        Theta(4) = Theta(4) + x(i) * sqrt_a_alpha * g(i);
        Theta(5) = Theta(5) + x(i) * b(i);
    end
    
    Theta(6) = L;  % Set the last value directly
end
