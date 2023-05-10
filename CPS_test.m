% Only for planar structure
% 
% 
% Upper structure
%
% thickness_upper [ Top , ... , ... , Closest to dipole ]
% ______________________                 -
%                                      z_plus
%         o       <--- Dipole            -
%                                      z_minus     
%                                           
% ______________________                 -
% 
% thickness_below [ Closest to dipole , ... , ... , Bottom ]
%
% Lower structure
%
%
% n: From top to bottom

j = 0;
for i = 1:length(up.thickness)
    if up.grating(i) == 1
        upt(i+j) = up.thickness(i)*up.duty(i);
        upt(i+j+1) = up.thickness(i)*(1-up.duty(i));
        j = j+1;
    else
        upt(i+j) = up.thickness(i);
    end
end

j = 0;
for i = 1:length(down.thickness)
    if down.grating(i) == 1
        downt(i+j) = down.thickness(i)*down.duty(i);
        downt(i+j+1) = down.thickness(i)*(1-down.duty(i));
        j = j+1;
    else
        downt(i+j) = down.thickness(i);
    end
end

thickness_upper = flip(upt(2:end));
thickness_below = downt(2:end);

n = [flip(up.n),down.n(2:end)];
lambda0 = param.wavelengths;
eta_emitter=1; 
du = 0.001;
u_upper = 2;

% LEE = 0;

z_plus = param.Org_thickness- param.Dipole_pos;
z_minus = param.Dipole_pos; 


u = 0:du:u_upper;
[K_TEh, K_TMh, F_out, theta_deg_air, eta_out, eta_wg_sub, eta_sp, eta_abs, eta_emitter_eff, eta_out_h, eta_wg_sub_h, eta_sp_h, eta_abs_h, eta_emitter_eff_h, eta_out_TMv, eta_wg_sub_TMv, eta_sp_TMv, eta_abs_TMv, eta_emitter_eff_TMv, u_crit, K_h, uK_h, K_h_prime, uK_h_prime, K_TMv, uK_TMv, K_TMv_prime, uK_TMv_prime, K, uK, K_prime, uK_prime, Purcell_h, Purcell_TMv, Purcell] = Func_CPS_clean(thickness_upper, thickness_below, z_plus, z_minus, n, lambda0, eta_emitter, du, u_upper);
