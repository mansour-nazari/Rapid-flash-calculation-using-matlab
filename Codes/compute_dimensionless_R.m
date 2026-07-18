function RD = compute_dimensionless_R(R, P, T, Rgas)
% COMPUTE_DIMENSIONLESS_R scales residual vector R to dimensionless form
%
% Inputs:
%   R     - residual vector [R1; R2; R3; R4; R5; R6] (6x1)
%   P     - pressure (scalar)
%   T     - temperature (scalar)
%   Rgas  - gas constant (scalar)
%
% Output:
%   RD    - dimensionless residual vector (6x1)

    scale_sqrt = sqrt(P) / (Rgas * T);
    scale_linear = P / (Rgas * T);

    RD = zeros(6,1);
    RD(1:4) = scale_sqrt * R(1:4);
    RD(5)   = scale_linear * R(5);
    RD(6)   = R(6);  % unscaled
end
