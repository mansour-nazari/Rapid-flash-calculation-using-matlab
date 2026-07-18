function [zl,zg]=Zfactor(E0G,E1G,E2G,E0L,E1L,E2L,Al,Bl,Ag,Bg,d1,d2)

coeff=[1 E2L E1L E0L];
Zl=roots(coeff);
zl=1;zg=1;           
num_real = 0; % we want to count the real roots
for number_root = 1: length (Zl)
    root = Zl (number_root);
    if root == real(root) % checking if the root is real or imaginary
        num_real = num_real + 1;
        real_roots(num_real) = root; %#ok     we save the real roots
    end
end 
    if num_real == 1 % if only one root is real, it is the answer (two conjugate roots)
        zl = real_roots (num_real);
    elseif num_real==3 % if three roots are real, we check the Gibs Energy
        % calling a function which gets two real roots and calculate Gibbs energy
        zl = Gibbs_energy (min(real_roots (:)), max(real_roots (:)), Al, Bl, d1, d2);
    end      
            
            % gas
            coeff1=[1 E2G E1G E0G];
            Zg=roots(coeff1);
            num_real = 0; % we want to count the real roots
for number_root = 1: length (Zg)
    root = Zg (number_root);
    if root == real(root) % checking if the root is real or imaginary
        num_real = num_real + 1;
        real_roots(num_real) = root; %ok  we save the real roots
    end
end 
    if num_real == 1 % if only one root is real, it is the answer (two conjugate roots)
        zg = real_roots (num_real);
    elseif num_real==3  % if three roots are real, we check the Gibs Energy
        % calling a function which gets two real roots and calculate Gibbs energy
        zg = Gibbs_energy (min(real_roots (:)), max(real_roots (:)), Ag, Bg, d1, d2);
    end   
            
            
            
end
            
 
