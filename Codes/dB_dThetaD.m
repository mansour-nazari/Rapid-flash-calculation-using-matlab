function dBdThetaD = dB_dThetaD(k)
% DBDTHETAD returns ∂B/∂Θ_Dk = δ_{5k}, i.e. 1 if k == 5, else 0
%
% Input:
%   k - index (1 to 6)
%
% Output:
%   dBdThetaD - scalar (0 or 1)

    dBdThetaD = double(k == 5);  % Kronecker delta δ_5k
end