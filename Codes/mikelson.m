clc
clear
close all

%% --- Compositions from Excel ---
[data, txt] = xlsread('pvt.xlsx','PROP'); % Component properties
CompList = txt(2:end,1);                  % Component names

%% --- inputs ---
P = 600;               % Pressure [psia]
T = 300 + 459.67;      % Temperature [°R]
z = data(:,1) / 100;   % Overall composition
Tc = data(:,3)+459.67; % Critical Temperatures [°R]
Pc = data(:,2);        % Critical Pressures [psia]
omega = data(:,4);     % Acentric factors
hi=data(:,end-1);      % h values
gi=data(:,end);         % g values
hi=zeros(length(hi),1);
gi=zeros(length(hi),1);
eos_flag = 'PR';       % 'PR' or 'SRK'
R = 10.73159; % Gas constant in ft^3·psi/(lbmol·°R)
N = length(z);

%%
% EOS selection
for i = 1:N
    Tr = T / (Tc(i));
    if strcmpi(eos_flag, 'PR')
        kappa = 0.37464 + 1.54226*omega(i) - 0.26992*omega(i)^2; % m(omega) 
        alpha(i) = (1 + kappa * (1 - sqrt(Tr)))^2; % alpha
        a(i) = 0.457235529 * R^2 * Tc(i)^2 / Pc(i); % a valves
        b(i) = 0.077796074 * R * Tc(i) / Pc(i); % b valvues
        d1 = 1 + sqrt(2);
        d2 = 1 - sqrt(2);
        u=2;
        w=-1;
    elseif strcmpi(eos_flag, 'SRK')
        kappa = 0.480 + 1.574*omega(i) - 0.176*omega(i)^2; % m(omega) 
        alpha(i) = (1 + kappa * (1 - sqrt(Tr)))^2; % alpha
        a(i) = 0.42748 * R^2 * Tc(i)^2 / Pc(i) ; % a valves
        b(i) = 0.08664 * R * Tc(i) / Pc(i); % b valvues
        d1 = 1;
        d2 = 0;
        u=1;
        w=0;
    else
        error('Unsupported EOS type');
    end
end
a=a'; alpha=alpha';b=b';

%%
K=wilsons(T,P,Tc-459.67,Pc,omega); % wilosn equ
nv=0.2; % nv initial guess
x = z./(1+nv*(K-1)); % x calc
x=x/sum(x); % normalize
y = K.*z./(1+nv*(K-1)); % y calc
y=y/sum(y); % normalize y

nL=1-nv; % liquid mol

Theta_L = compute_theta_L(a, alpha, hi, gi, x, b, nL); % theta liquid equations 
Theta_F = compute_theta_L(a, alpha, hi, gi, z, b, 1); % theta feed
 
tol=1;
index=0;
success = false;

%% main loop
while tol>10^-13
    index=1+index;
    Theta_V=(Theta_F-nL*Theta_L)/(1-nL); % vapor theta
    Theta_V(end)=1-Theta_L(end);
    
    [Psi_L, a_m_L, b_m_L] = compute_psi_am_bm_from_theta(a, alpha, hi, gi, Theta_L); % liquid psi,a_mixture, b_mixture
    [Psi_V, a_m_V, b_m_V] = compute_psi_am_bm_from_theta(a, alpha, hi, gi, Theta_V); % vapor psi,a_mixture, b_mixture

    % EOS A,B for gas and liquid
    Al = a_m_L*P/(R*T)^2;
    Bl = (b_m_L)*P/(R*T);
    Ag = a_m_V*P/(R*T)^2;
    Bg = (b_m_V)*P/(R*T);

    % EOS coefficeint
    [E0G,E1G,E2G,E0L,E1L,E2L]=EOS_coeff(Ag,Bg,Al,Bl,u,w);
    % EOS ROOT SELECTION USING GIBBS ENERGY EQU FOR LIQUID AND VAPOR
    [zl,zg]=Zfactor(E0G,E1G,E2G,E0L,E1L,E2L,Al,Bl,Ag,Bg,d2,d1);
    
    % Bi and Aij for fugacity calculation
    Bi=P/(R*T)*b;
    Aij_L=P/(R*T)^2*Psi_L;
    Aij_V=P/(R*T)^2*Psi_V;
    
    % calculate fugacity for vapor and liquid
    ln_phi_L = compute_ln_phi_L_vector(zl, Al, Bl, Bi, Aij_L, d1, d2);
    ln_phi_V = compute_ln_phi_L_vector(zg, Ag, Bg, Bi, Aij_V, d1, d2);
    phi_l=exp(ln_phi_L);
    phi_v=exp(ln_phi_V);
    
    % updating k value
    K_n=phi_l./phi_v;
    K=K_n;
    x = z./(1+Theta_V(end)*(K-1)); % x calc
    x=x/sum(x);
    y = K.*z./(1+Theta_V(end)*(K-1)); % y calc
    y=y/sum(y);
    % 6 R equation calc
    R1 = compute_R_vector(a, alpha, hi, gi, b, x, Theta_L, z, K);
    
    % Ai and Bi and dA/dtheta
    [Ai, Bi, dAij_dTheta] = Jac_parametrs(a, alpha, b, hi, gi, P, T, R);
    
    % dimensionless theta for vapor and liquid
    Theta_D_L = compute_dimensionless_Theta(Theta_L, P, T, R);
    Theta_D_V = compute_dimensionless_Theta(Theta_V, P, T, R);

    % dA/dtheta
    dA_dThetaD_L = dA_dThetaD_vector(Theta_D_L);
    dA_dThetaD_V = dA_dThetaD_vector(Theta_D_V);
    
    %dB/dtheta
    dB_dThetaD_L = zeros(1,6); dB_dThetaD_L(5) = 1;
    dB_dThetaD_V = zeros(1,6); dB_dThetaD_V(5) = 1;
    
    %dz/dtheta
    dZ_dThetaD_L = compute_dZ_dThetaD(zl, Al, Bl, d1, d2, dA_dThetaD_L, dB_dThetaD_L);
    dZ_dThetaD_V = compute_dZ_dThetaD(zg, Ag, Bg, d1, d2, dA_dThetaD_V, dB_dThetaD_V);
    
    % dln(phi)/dtjeta
    dlnphi_L = dlnphi_dThetaD(zl, Al, Bl, Bi, Aij_L, dZ_dThetaD_L, dA_dThetaD_L, dB_dThetaD_L, dAij_dTheta, d1, d2);
    dlnphi_V = dlnphi_dThetaD(zg, Ag, Bg, Bi, Aij_V, dZ_dThetaD_V, dA_dThetaD_V, dB_dThetaD_V, dAij_dTheta, d1, d2);

    % dtheta(V)/dtheta(L)
    dThetaV_dThetaL = dThetaDv_dThetaDl_matrix(Theta_D_L,Theta_D_V);

    % dk/dtheta
    dK_dThetaDl = compute_dK_dThetaDl_matrix(K, dlnphi_L, dThetaV_dThetaL);

    % dx/dtheta and dy/dtheta
    [dx_dTheta, dy_dTheta] = compute_dx_dy_dTheta(K, z, dK_dThetaDl, Theta_D_L(end));

    % jacobian matrix
    Jac = compute_main_jacobian(Ai, hi, gi, Bi, dx_dTheta, dy_dTheta);

    % dimensionless R calc
    RD = compute_dimensionless_R(R1, P, T, R);

    % updating theta using newton equ
    Theta_D_new = (Theta_D_L' - Jac\RD)';
    
    % convert dimensionless theta to theta
    Theta_L = convert_ThetaD_to_Theta(Theta_D_new, P, T, R);

    nL=Theta_L(end);
    tol=norm(R1);
    
        % Convergence check
    if tol < 10^-6
        success = true;
        break;
    end
end


%% --- Results ---
if success
    disp('Rapid Flash converged:');
    disp(eos_flag);
    fprintf('P = %.6f\n\n', P);
    fprintf('T = %.6f\n\n', T);
    fprintf('nL = %.6f\n\n', nL);
    fprintf('Component\t\tx\t\ty\t\tK\n');
    for i = 1:N
        fprintf('%-15s\t%.5f\t%.5f\t%.5f\n', CompList{i}, x(i), y(i),K(i));
    end
else
    disp('Flash did not converge.');
end
