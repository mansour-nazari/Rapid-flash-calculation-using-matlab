function dK_dThetaDl = compute_dK_dThetaDl_matrix(K, dlnphi_V, dThetaV_dThetaL)
% COMPUTE_DK_DTHETADL_MATRIX computes ∂K_i/∂Θ_Dk^L for all i and k
%
% Inputs:
%   K                      - N×1 vector of K_i values
%   dlnphiV_dThetaDv_all   - N×6 matrix: row i is ∂ln(φ_i^V)/∂Θ_Dl^V (1×6)
%   dThetaDv_dThetaDl      - 6×6 Jacobian matrix: ∂Θ_Dl^V / ∂Θ_Dk^L
%   dlnphiV_dThetaDl_all   - N×6 matrix: row i is ∂ln(φ_i^V)/∂Θ_Dk^L
%
% Output:
%   dK_dThetaDl            - N×6 matrix of ∂K_i/∂Θ_Dk^L
    N = length(K);
    dK_dThetaDl = zeros(N, 6);

    for i = 1:N
        for k = 1:6
            term=0;
            for l=1:6
                term=dlnphi_V(i,l)*dThetaV_dThetaL(l,k)-dlnphi_V(i,k)+term;
            end
            dK_dThetaDl(i,k)=K(i)*term;
        end
    end
end
