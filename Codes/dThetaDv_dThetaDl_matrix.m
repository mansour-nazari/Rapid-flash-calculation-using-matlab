function dThetaDv_dThetaDl = dThetaDv_dThetaDl_matrix(ThetaDl, ThetaDv)
% DTHETADV_DTHETADL_MATRIX returns the 6x6 Jacobian matrix of
% ∂Θ_Dl^V / ∂Θ_Dk^L as defined by the piecewise equation.
%
% Inputs:
%   ThetaDl - 6x1 vector of Θ_D^L values
%   ThetaDv - 6x1 vector of Θ_D^V values
%
% Output:
%   dThetaDv_dThetaDl - 6x6 matrix of derivatives

    dThetaDv_dThetaDl = zeros(6, 6);
    TL6 = ThetaDl(6);

    for l = 1:6
        for k = 1:6
            if l == k && k <= 5
                dThetaDv_dThetaDl(l, k) = -(1 - TL6) / TL6;
            elseif l == 6 && k <= 5
                dThetaDv_dThetaDl(l, k) = (ThetaDl(k) - ThetaDv(k)) / (1 - ThetaDl(k))^2;
            elseif l == 6 && k == 6
                dThetaDv_dThetaDl(l, k) = -1;
            else
                dThetaDv_dThetaDl(l, k) = 0;
            end
        end
    end
end
