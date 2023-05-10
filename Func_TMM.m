%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4th floor                         %
% interface (count = 3)             %
% 3rd floor                         %
% interface (count = 2)             %
% 2nd floor                         %
% interface (count = 1)             %
% 1st floor                         %
%                                   %
% Direction of light : Up           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TMM([0 75e-9 50e-9 0], [1.74 2 3 1], 500*1e-9, 0.5)

% r : reflectance
% t : transmittance
% R : energy reflectivity
% T : energy transmittance

function [r_fin_TE,R_fin_TE,t_fin_TE,T_fin_TE,r_fin_TM,R_fin_TM,t_fin_TM,T_fin_TM] = Func_TMM(thickness, n, lambda0, u) % Vector / Vector / scalar / scalar
floor = length(thickness);
n_init = n(1);
n_fin = n(length(n));

k = n.*2*pi/lambda0;
k0 = 2*pi/lambda0;

%%TE
    M_0N_TE = eye(2);
    for count = 1:floor-1
        
        r_TE(count) = ( n(count)*sqrt(1-n_init^2/n(count)^2*u^2) - n(count+1)*sqrt(1-n_init^2/n(count+1)^2*u^2) ) / ( n(count)*sqrt(1-n_init^2/n(count)^2*u^2) + n(count+1)*sqrt(1-n_init^2/n(count+1)^2*u^2) );
        t_TE(count) = 1 + r_TE(count);
        
        Mij_TE_11(count) = (1/t_TE(count));
        Mij_TE_12(count) = (1/t_TE(count))*r_TE(count);
        Mij_TE_21(count) = (1/t_TE(count))*r_TE(count);
        Mij_TE_22(count) = (1/t_TE(count));
        
        Mi_TE_11(count) = exp(-i*k0*n_init*sqrt(n(count)^2/n_init^2-u.^2)*thickness(count));%
        Mi_TE_12(count) = 0;
        Mi_TE_21(count) = 0;
        Mi_TE_22(count) = exp(+i*k0*n_init*sqrt(n(count)^2/n_init^2-u.^2)*thickness(count));%
        
        Mij_TE(:,:,count) = [Mij_TE_11(count) Mij_TE_12(count) ; Mij_TE_21(count) Mij_TE_22(count)];
        Mi_TE(:,:,count) = [Mi_TE_11(count) Mi_TE_12(count) ; Mi_TE_21(count) Mi_TE_22(count)];

    end

    for count = 1 : floor-1
        if count == floor-1
            M_0N_TE =  M_0N_TE * Mij_TE(:,:,count);
        else
            M_0N_TE = M_0N_TE * Mij_TE(:,:,count) * Mi_TE(:,:,count+1);
        end
    end
    
    r_fin_TE = M_0N_TE(2,1)/M_0N_TE(1,1);
    t_fin_TE = 1/M_0N_TE(1,1);
    R_fin_TE = (abs(r_fin_TE))^2;
    T_fin_TE = (abs(t_fin_TE))^2 * n_fin/n_init * sqrt(1-n_init^2/n_fin^2*u.^2) / sqrt(1-u.^2);
    
%%TM
    M_0N_TM = eye(2);
    for count = 1:floor-1
        
        r_TM(count) = - ( n(count)*sqrt(1-n_init^2/n(count+1)^2*u^2) - n(count+1)*sqrt(1-n_init^2/n(count)^2*u^2) ) / ( n(count)*sqrt(1-n_init^2/n(count+1)^2*u^2) + n(count+1)*sqrt(1-n_init^2/n(count)^2*u^2) ); % - : Centrioni
        t_TM(count) = n(count)/n(count+1) * (1+r_TM(count));

        Mij_TM_11(count) = (1/t_TM(count));
        Mij_TM_12(count) = (1/t_TM(count))*r_TM(count);
        Mij_TM_21(count) = (1/t_TM(count))*r_TM(count);
        Mij_TM_22(count) = (1/t_TM(count));
        
        Mi_TM_11(:,count) = exp(-i*k0*n_init*sqrt(n(count)^2/n_init^2-u.^2)*thickness(count));
        Mi_TM_12(:,count) = 0;
        Mi_TM_21(:,count) = 0;
        Mi_TM_22(:,count) = exp(+i*k0*n_init*sqrt(n(count)^2/n_init^2-u.^2)*thickness(count));
        
        Mij_TM(:,:,count) = [Mij_TM_11(count) Mij_TM_12(count) ; Mij_TM_21(count) Mij_TM_22(count)];
        Mi_TM(:,:,count) = [Mi_TM_11(count) Mi_TM_12(count) ; Mi_TM_21(count) Mi_TM_22(count)];

    end

    for count = 1 : floor-1
        if count == floor-1
            M_0N_TM =  M_0N_TM * Mij_TM(:,:,count);
        else
            M_0N_TM = M_0N_TM * Mij_TM(:,:,count) * Mi_TM(:,:,count+1);
        end
    end
    
    r_fin_TM = M_0N_TM(2,1)/M_0N_TM(1,1);
    t_fin_TM = 1/M_0N_TM(1,1);
    R_fin_TM = (abs(r_fin_TM))^2;
    T_fin_TM = (abs(t_fin_TM))^2 * n_fin/n_init * sqrt(1-n_init^2/n_fin^2*u.^2) / sqrt(1-u.^2);
end