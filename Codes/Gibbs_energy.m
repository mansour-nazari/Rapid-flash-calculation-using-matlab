function z = Gibbs_energy (Z_low, Z_high, A, B, delta1, delta2)
% This function calculate the real root of EOS equation.
% The root which have lower Gibbs energy is favorable.
% Input (lower root, higher root, A, B, delta1, delta2)
% Attention: the middle root is incorrect so we dont use it (dp/dv<0)
% Output: compressibility factor
DGRT = (Z_high-Z_low) - log ((Z_high - B)/(Z_low - B)) - (A/ (B*(delta2 - delta1)))*log (((Z_high + delta2*B)*(Z_low + delta1*B))/((Z_high + delta1*B)*(Z_low + delta2*B)));
if DGRT < 0
    z = Z_high;
else
    z = Z_low;
end
end
