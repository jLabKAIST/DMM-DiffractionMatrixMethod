function [rd,td,struc] = TMM_DMM2D(param,up)

% RCWA_DOWN PARAMETERS

uxsize = param.uxsize;
uysize = param.uysize;
wavelengths = param.wavelengths;

struc = zeros(round(sum(up.thickness*1e9)),100,100);
numlayer = length(up.grating);

rd.TETE = zeros(length(uxsize), length(uysize),13,13);
rd.TETM = zeros(length(uxsize), length(uysize),13,13);
rd.TMTE = zeros(length(uxsize), length(uysize),13,13);
rd.TMTM = zeros(length(uxsize), length(uysize),13,13);

rTETE = zeros(length(uxsize), length(uysize));
rTMTM = zeros(length(uxsize), length(uysize));

td.TETE = zeros(length(uxsize), length(uysize),13,13);
td.TETM = zeros(length(uxsize), length(uysize),13,13);
td.TMTE = zeros(length(uxsize), length(uysize),13,13);
td.TMTM = zeros(length(uxsize), length(uysize),13,13);

tTETE = zeros(length(uxsize), length(uysize));
tTMTM = zeros(length(uxsize), length(uysize));

thickness_prime = up.thickness;
n_prime = up.n;

k = 0; 
for i = 1:numlayer
    struc(k+1:round(k+up.thickness(i)*1e9),:,:) = up.n(i);
    k = k+round(up.thickness(i)*1e9);
end


for uxindex = 1:length(uxsize)
parfor uyindex = 1:length(uysize)
     
ux = uxsize(uxindex);
uy = uysize(uyindex);

u_temp = sqrt(ux^2 + uy^2);

[r_TE_temp,R_TE_temp,t_TE_temp,T_TE_temp,r_TM_temp,R_TM_temp,t_TM_temp,T_TM_temp] = Func_TMM(thickness_prime, n_prime, wavelengths, u_temp);
    
rTETE(uxindex,uyindex) = r_TE_temp;
rTMTM(uxindex,uyindex) = r_TM_temp;

tTETE(uxindex,uyindex) = t_TE_temp;
tTMTM(uxindex,uyindex) = t_TM_temp;

% Tup_amp_TETE0(uxindex,uyindex) = t_TE_temp;
% Tup_amp_TMTM0(uxindex,uyindex) = t_TM_temp;


retio;
end
end

rd.TETE(:,:,7,7) = rTETE;
rd.TMTM(:,:,7,7) = rTMTM;

td.TETE(:,:,7,7) = tTETE;
td.TMTM(:,:,7,7) = tTMTM;

end
