function [E0G,E1G,E2G,E0L,E1L,E2L]=EOS_coeff(Ag,Bg,Al,Bl,u,w)

            % gas phase
                E0G=-Ag*Bg-w*Bg^2-w*Bg^3;
                E1G=Ag+(w-u)*Bg^2-(u)*Bg;
                E2G=(u-1)*Bg-1;
            % liquid phase
                E0L=-Al*Bl-w*Bl^2-w*Bl^3;
                E1L=Al+(w-u)*Bl^2-(u)*Bl;
                E2L=(u-1)*Bl-1;
                
end