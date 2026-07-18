function Theta_D = compute_dimensionless_Theta(Theta, P, T, Rgas)
% COMPUTE_DIMENSIONLESS_THETA scales Theta vector to dimensionless form
%
% Inputs:
%   Theta - 6x1 vector of Θ_L values
%   P     - pressure (scalar)
%   T     - temperature (scalar)
%   Rgas  - gas constant (scalar)
%
% Output:
%   Theta_D - dimensionless Θ vector (6x1)

    scale_sqrt = sqrt(P) / (Rgas * T);
    scale_linear = P / (Rgas * T);

    Theta_D = zeros(1,6);
    Theta_D(1:4) = scale_sqrt * Theta(1:4);
    Theta_D(5)   = scale_linear * Theta(5);
    Theta_D(6)   = Theta(6);  % unscaled
end
