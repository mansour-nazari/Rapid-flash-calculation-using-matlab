function Theta = convert_ThetaD_to_Theta(ThetaD, P, T, R)
% CONVERT_THETAD_TO_THETA rescales dimensionless Θ_D back to dimensional Θ
%
% Inputs:
%   ThetaD - 6x1 vector of dimensionless Theta_D
%   P      - pressure (scalar)
%   T      - temperature (scalar)
%   R      - gas constant (scalar)
%
% Output:
%   Theta  - 6x1 vector of dimensional Theta

    scale_sqrt = sqrt(P) / (R * T);
    scale_linear = P / (R * T);

    Theta = zeros(1,6);
    Theta(1:4) = ThetaD(1:4) / scale_sqrt;
    Theta(5)   = ThetaD(5) / scale_linear;
    Theta(6)   = ThetaD(6);  % unscaled
end
